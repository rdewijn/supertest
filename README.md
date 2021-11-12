# Supertest

##### Description

The `Supertest` shiny operator is for mapping peptides to upstream kinases and performing a test on differential kinase activity based on peptide sets.

##### Usage

Input projection|.
---|---
`x-axis`        | type, description 
`y-axis`        | type, description 
`row`           | type, description 
`column`        | type, description 
`colors`        | type, description 
`labels`        | type, description 

Output relations|.
---|---
`Operator view`        | view of the Shiny application

##### Details

The R-package [globaltest](https://www.bioconductor.org/packages/release/bioc/html/globaltest.htmlpackage) is aplied to perform a hypothesis test on a set of peptides.

Goeman JJ, Oosting J (2020). Globaltest: testing association of a group of genes with a clinical variable. R package version 5.48.0.

Goeman JJ, van de Geer SA, van Houwelingen JC (2006). “Testing against a high-dimensional alternative.” Journal of the Royal Statistical Society, Series B, 477-493.

Goeman JJ, van de Geer SA, de Kort F, van Houwelingen JC (2004). “A global test for groups of genes: testing association with a clinical outcome.” Bioinformatics, 93-99.

##### See Also


