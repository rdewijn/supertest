library(shiny)
library(tercen)
library(tidyverse)
library(globaltest)
library(reshape2)
library(data.table)
library(ggsci)
library(pgSupertest)

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
  
  mode = reactive({
    getMode(session)
  })
  
  
  
  observe({
    dfin = dataIn()
    
    if(!isRunView(mode()) ){
      shinyjs::disable("done")
    }
    
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
        ungroup() %>%
        mutate(fdr = p.adjust(p, "fdr"))%>%
        select(Kinase_Name, N, delta, p, fdr, Family)
      removeNotification(id)
      return(result)
    })
    
    getGrp = reactive({
      grp = dfin$grp %>%
        as.factor()
    })
    
    grpText = reactive({
      grp = getGrp()
       txt = paste("Upstream kinase test using a grouping factor with levels", levels(grp)[1], "and", levels(grp)[2])
    })
    
    upText = reactive({
      grp = getGrp()
      txt = paste("Higher activity in the", levels(grp)[2],"group")
    })
    
    downText = reactive({
      grp = getGrp()
      txt = paste("Lower activity in the", levels(grp)[2],"group")
    })
    
    output$kup = renderTable({
        getResultTable() %>% 
        filter(delta >= 0)
    })
    
    output$kdn = renderTable({
        getResultTable() %>%
        filter(delta < 0)
    })
    
    output$grpTxt = renderText({
      grpText()
    })
    
    output$uptext = renderText({
      upText()
    })
    
    output$downtext = renderText({
      downText()
    })
    
    observeEvent(input$done, {
      ctx <- context()
      getResultTable() %>%
        ungroup() %>%
        mutate(.ri = 0:(n()-1), .ci = 0) %>%
        select(.ri, .ci ,p) %>%
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
    stop("Need exactly one color for the grouping")
  }
  if (!("ID" %in% ctx$rnames[[1]]) )stop("Factor ID for identifying spots is required")
  id = data.frame(.ri =unique(df$.ri), ID=ctx%>%rselect("ID"))
  
  df = df %>% 
    left_join(id, by = ".ri") %>%
    mutate( across(!where(is.numeric), as.factor)) %>%
    arrange(.ri, .ci)
}

getMode = function(session){
  # retreive url query parameters provided by tercen
  query = parseQueryString(session$clientData$url_search)
  return(query[["mode"]])
}

isRunView <- function(mode) {
  !is.null(mode) && mode == "run"
}

