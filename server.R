library(shiny)
library(tercen)
library(tidyverse)
library(globaltest)
library(reshape2)
library(data.table)
library(ggsci)

############################################
#### This part should not be modified
getCtx <- function(session) {
  # retreive url query parameters provided by tercen
  query <- parseQueryString(session$clientData$url_search)
  token <- query[["token"]]
  taskId <- query[["taskId"]]
  
  # create a Tercen context object using the token
  ctx <- tercenCtx(taskId = taskId, authToken = token)
  return(ctx)
}
####
############################################

server <- shinyServer(function(input, output, session) {
  
  db = reactive({
    upstreamDb()
  })
  
  s1 = reactive({
    TableS1()
  })
  
  dataIn = reactive({
    getValues(session)
  })
  
  
  
  
  observe({
    dfin = dataIn()
    
    selectMappings = reactive({
      tabs1 = s1() %>%
        filter(GroupName == "TK")
      db() %>%
        ungroup() %>% 
        filter(PepProtein_SeqHomology >= input$seqhom) %>%
        filter(family == "PTK") %>% 
        filter(Kinase_Name %in% tabs1$KR_Name) %>%
        filter(ID %in% dfin$ID)
    })
    
    scores = reactive({
      dbi = selectMappings() %>% 
        scoreIviv(input$iviv)
      dbp = selectMappings() %>%
        scorePNet(ranks = input$pn.breakpoints, scores = c(input$pn.low, input$pn.low, input$pn.high))
      combinedScores(list(dbi, dbp), input$p0) %>%
        normalizeScores()
    })
    
    kinasetest = reactive({
      X = acast(dfin, .ci ~ ID, value.var = ".y")
      grp = acast(dfin, .ci ~ ID, value.var = "grp")[,1]
      result = scores() %>%
        filter(sc.nor > input$thr) %>%
        group_by(Kinase_Name) %>%
        do(gtest(.,X, grp))
    })
    
    getResultTable = reactive({
      id = showNotification(paste("Please wait ..."), duration = NULL, closeButton = FALSE, type = "message")
      result =
        kinasetest() %>%
        filter(N >= input$pss) %>%
        arrange(p) %>%
        left_join(s1(), by = c(Kinase_Name = "KR_Name")) %>%
        mutate(fdr = p.adjust(p, "fdr"))%>%
        select(Kinase_Name, N, delta, p, fdr, Family)
      removeNotification(id)
      return(result)
    })
    
    output$kup = renderTable({
        getResultTable() %>% 
        filter(delta >= 0)
    })
    
    output$kdn = renderTable({
        getResultTable() %>%
        filter(delta < 0)
    })
    
    output$vplot = renderPlot({
      result = getResultTable()  
      p = ggplot(result, aes(x = delta, y = -log10(fdr), label = Kinase_Name, size = N, colour =Family ))
        p = p + geom_text() + theme_bw()
        p = p + xlab("delta") + ylab("-log10(fdr)") + guides(color  = FALSE)
        #p = p + xlim( max(abs(result$delta)) *c(-1.1, 1,1) ) + ylim( c(0, max(-log10(result$fdr))) )
        print(p)
    })
    
    observeEvent(input$done, {
      ctx <- context()
      getResultTable() %>%
        ungroup() %>%
        mutate(.ri = 0:(n()-1), .ci = 0) %>%
        ctx$addNamespace() %>%
        ctx$save()
    })
    context <- reactive({
      getCtx(session)
    })
  })
  
})


getValues <- function(session){
  ctx <- getCtx(session)
  df <- ctx %>% 
    select(.y, .ri, .ci)
  if(length(ctx$colors) == 1){
    df = df %>% bind_cols(data.frame(grp = ctx$select(ctx$colors)[[1]]))
  } else {
    stop("Need axctly one color for the grouping")
  }
  if (!("ID" %in% ctx$rnames[[1]]) )stop("Factor ID for identifying spots is required")
  id = data.frame(.ri =unique(df$.ri), ID=ctx%>%rselect("ID"))
  
  df = df %>% 
    left_join(id, by = ".ri") %>%
    mutate( across(!where(is.numeric), as.factor)) %>%
    arrange(.ri, .ci)
}

upstreamDb = function(){
  load("./d/180509_86312_86402_87102_UpstreamDb.RData")
  UpstreamDatabase
}

TableS1 = function(){
  df = read.delim("./d/TableS1 Kinases in Kinome Render.txt")
}

scoreIviv = function(db, score){
  db %>% 
    filter(Database == "iviv") %>%
    mutate(s = score)
}

scorePNet = function(db, ranks, scores){
  db %>% 
    mutate(
      s  = case_when(
        Kinase_Rank > ranks[2] ~ scores[3],
        Kinase_Rank > ranks[1] ~ scores[2],
        TRUE ~ scores[1])
    )
}

combinedScores = function(dblist, minscore){
  dbc = bind_rows(dblist) %>%
    group_by(Kinase_Name, ID) %>%
    summarise(sc = 1-prod(1-s))
  
  db.all = expand.grid(levels(droplevels(dbc$Kinase_Name)), levels(droplevels(dbc$ID)) )
  colnames(db.all) = c("Kinase_Name", "ID")
  
  db.all %>%
    left_join(dbc, by = c("Kinase_Name", "ID")) %>%
    mutate(sc.final = case_when(
      is.na(sc) ~ minscore,
      sc < minscore ~ minscore,
      TRUE ~ sc))
}

normalizeScores = function(db){
  sumdf = db %>% group_by(ID) %>% summarise(sumsc = sum(sc.final))
  db %>%
    left_join(sumdf, by = "ID") %>%
    group_by(ID) %>%
    mutate(sc.nor = sc.final/sumsc)
}

gtest = function(kinase.set, X, grp){
  bset = colnames(X) %in% kinase.set$ID
  X = X[, bset,drop = FALSE]
  if(dim(X)[2]>0){
    globtest = gt(grp ~ X, directional = TRUE, standardize = TRUE)
    p = attr(globtest, "result")[1,1]
    delta = mean(gdelta(X, grp))
  } else {
    p = NaN
    delta = NaN
  }
  result = data.table(p = p, N = dim(X)[2], delta = delta, gt = list(globtest))
}

gdelta = function(Xin, grp){
  X1 = Xin[grp == levels(grp)[1],,drop = FALSE]
  X2 = Xin[grp == levels(grp)[2],,drop= FALSE]
  M1 =apply(X1, 2, mean)
  M2 =apply(X2, 2, mean)
  delta = M2 - M1
}

