# QC plots
plot_hist_perc <- function ( x, plot_this="percent.mt", ccc="#6e9be5" , cut_here = 15, msg = "Cells with MT >", binwidth= 5 , y=.5) {
  
   
 x$x1 = x[[plot_this]]
  
 ggg= ggplot(x, aes(x=x1)) + 
    geom_histogram(aes(y=(..count..)/sum(..count..)), colour='black' , fill=ccc, alpha=.2 , binwidth= binwidth )+
    # geom_density(alpha=0, fill="grey") +
    theme(legend.position="none",  legend.key = element_blank(),
          # element_blank()
          axis.text.y = element_text(size= 25 ),
          axis.text.x =  element_text(size= 25 ),
          axis.title.x = element_text(size=25),
          axis.title.y     = element_text(size=25), 
          legend.text      =element_text(size=25),
          legend.title = element_text(size=25),
          plot.title = element_text(size = 40, face = "bold", hjust = 0.5)
          # hjust centers the title
    ) + theme(panel.grid.major = element_blank()
              , panel.grid.minor = element_blank()
              ,panel.background = element_blank()
              , axis.line = element_line(colour = "black") # plot border
    ) + ylab ( "percent") + 
   geom_vline(xintercept = cut_here, linetype = 2, color = "red") +
   annotate(geom = "text", x = cut_here , y = y, label = paste0 ( msg , cut_here ), color = "red",
            angle = 90, vjust=2.5) + xlab (plot_this)
  
  
  return ( ggg )
  
  
}

# pahway 
get.enrich4 <- function (g,
                         refdb 
                         ,g.name="GENE_SYMBOL" # column for input dataframe for gene name 
                         , m.name="external_gene_name" # reference column for gene
                         # the two above will be merged so that we can get the entrezgene id 
                         , e.name="external_gene_name"
                         , logFC = "logFC" # name of the column to subset up and down
                         , species="Hs"
                         ,bg # usually all genes 
                         , fdr=.05
                         , goana_db = NA 
                         , triml = 60 # used to trim the names ince it can get really freq long. 
                         , custom = 'KEGG' 
                         , type = c("KEGG", "HALLMARK", "REACT","BIOC", "BP","MF","CC")
                         , total_annt = 0 # 
){
  
  # make a pseudo translator 
  decode = goana_db
  decode$GeneID = decode$PathwayID
  g <- merge(g,refdb, by.x=g.name, by.y=m.name )
  
  g = g[!duplicated(g[[g.name]]), ]
  
  g = g[!is.na ( g$entrezgene_id ) , ]
  
  # this function is used to trim the names  
  tm.this <- function (x){
    x$Term <- strtrim(x$Term,60)
    return (x)
  }
  
  up <- unique ( g[g[[logFC]] > 0 ,e.name] )
  down <-unique (  g[g[[logFC]] < 0 ,e.name] )
  
  # bg is now force to have input so you need to provide a background! 
  bg = merge(data.frame(bg=bg),refdb, by.x="bg", by.y=m.name )
  bg <- unique ( refdb[ , e.name] )
  
  bg = bg[ ! is.na( bg )]
  
  
  # run hypergeometric test 
  
  kegg <- kegga(list(Up=up, Down=down)
                , gene.pathway=goana_db
                , pathway.names = decode
                ,species=species,
                universe=bg, FDR=fdr
                
  )
  
  raw = kegg 
  # make sure the fdr cutoff is correct 
  kegg = kegg[(kegg$P.Up < .05 | kegg$P.Down < .05 ), ]
  
  
  
  keggUP = kegg[order(kegg$P.Up), ]
  keggDown = kegg[order(kegg$P.Down), ]
  kegg$p = ifelse ( kegg$P.Up < kegg$P.Down, kegg$P.Up, kegg$P.Down)
  # the last one makes sure that p value is not specific to up or down rather it capture both and you can subset later
  keggn = kegg[ order(kegg$p), ]
  keggn$P.Up = ifelse ( keggn$P.Up < .05, "*", "ns")
  keggn$P.Down = ifelse ( keggn$P.Down < .05, "*", "ns")
  
  keggn = keggn[keggn$p < .05, ] 
  
  
  # figure genes for enrichments 
  
  ingene = data.frame ( stringsAsFactors = F )
  
  # p="BP_BIOLOGICAL_ADHESION"
  # get genes and ratios associated with each pathways 
  for ( p in keggn$Pathway) {
    
    
    ga = unique ( goana_db[goana_db$PathwayID == p, ]$GeneID )
    inset = intersect ( unique ( c(up) ), ga ) 
    ratio = round ( length ( inset ) / length ( ga ) , 2)
    GeneRatio = paste0( length ( inset )  , "/", length ( ga )  )
    geneID = inset 
    geneName = ( refdb[ refdb[[ e.name ]] %in% inset,  m.name ]  )
    geneName = geneName[ !duplicated ( geneName )]
    
    ugene = geneName 
    
    inset = intersect ( unique ( c(down) ), ga ) 
    ratio = round ( length ( inset ) / length ( ga ) , 2)
    GeneRatio = paste0( length ( inset )  , "/", length ( ga )  )
    geneID = inset 
    geneName = ( refdb[ refdb[[ e.name ]] %in% inset,  m.name ]  )
    geneName = geneName[ !duplicated ( geneName )]
    
    dgene = geneName 
    
    
    
    
    
    
    ingene = rbind ( ingene, data.frame ( Pathway = p, GeneRatio = ratio
                                          , geneID=  toString(  paste ( geneID, collapse="," ) ) 
                                          , Up_regulated= toString(  ugene )
                                          , Down_regulated= toString(dgene ), stringsAsFactors = F
    ), stringsAsFactors = F 
    )
    
  }
  
  keggn = merge ( keggn, ingene, by="Pathway")
  row.names ( keggn ) = keggn$Pathway
  
  ########## end gene analysis 
  
  
  # you are allow at least 1 custom 
  
  type = unique ( c  ( custom, type ) )
  type = type [ type != "" ]
  
  enrich = list ()
  plots = list()
  
  for ( e in type ){
    enrich[[e]] = keggn[ grepl(paste0 ( "^", e   ), row.names ( keggn) ), ]
    if ( nrow ( enrich[[e]] ) > 0 ){
      
      enrich[[e]]$GeneRatio_Down =  round ( enrich[[e]]$Down / enrich[[e]]$N, 2) 
      enrich[[e]]$GeneRatio_Up =  round ( enrich[[e]]$Up / enrich[[e]]$N, 2) 
      # prioritize up always 
      enrich[[e]]$GeneRatio =  ifelse ( enrich[[e]]$P.Up == "*", enrich[[e]]$GeneRatio_Up, enrich[[e]]$GeneRatio_Down  )
      enrich[[e]]  = enrich[[e]] [ order ( enrich[[e]] $p, -enrich[[e]]$GeneRatio) , ]
      if ( total_annt > 0 ){
        plots[[e]] = plot.hyper ( enrich[[e]] , title = e, goana_db = goana_db , up=up, down=down, refdb=refdb, total=total_annt ) 
      }else {
        plots[[e]] = plot.hyper ( enrich[[e]] , title = e)
      }
    }
  }
  
  
  return ( list ( df=enrich
                  ,plots=plots  ) )
}


plot.hyper = function ( dfg, thres = 8 , title = "", path_length=50, total = 25, goana_db=NA, up=NA, down=NA , refdb=NA ){
  
  dfg = dfg[ order ( dfg$p, -  dfg$GeneRatio ), ]
  
  
  dfg.up  = dfg [ dfg$P.Up == "*", ]
  dfg.up$GeneRatio = dfg.up$Up/dfg.up$N * 100 
  dfg.up = dfg.up[ order ( dfg.up$p, -dfg.up$GeneRatio ), ]
  if ( nrow ( dfg.up) > 0 ){
    dfg.up$group = "up"
  }else{
    dfg.up = cbind ( dfg.up, data.frame ( group=as.character())   )
  }
  
  
  dfg.down  = dfg [ dfg$P.Down == "*", ]
  
  if ( nrow ( dfg.down) > 0 ){
    
    
    dfg.down$GeneRatio = dfg.down$Down/dfg.down$N * 100 
    dfg.down = dfg.down[ order ( dfg.down$p, -dfg.down$GeneRatio ), ]
    dfg.down$group = "down"
  }else{
    dfg.down = cbind ( dfg.down, data.frame ( group=as.character())   )
  }
  
  
  
  # get top up and down 
  total_down = ifelse(nrow ( dfg.down) <= total, total, nrow ( dfg.down) )
  total_up = ifelse(nrow ( dfg.up) <= total, total, nrow ( dfg.up) )
  
  
  dfg_full = rbind (  dfg.up[ , c("Pathway","p","GeneRatio", "group")], total_up 
                      , dfg.down[ , c("Pathway","p","GeneRatio", "group")], total_down) 
  
  
  dfg_full = dfg_full[ order (dfg_full$group , -dfg_full$p, abs ( dfg_full$GeneRatio), decreasing = T), ]
  dfg_full$Pathway = str_trunc ( as.character ( dfg_full$Pathway) , path_length )
  dfg_full$Pathway = make.unique(dfg_full$Pathway )
  
  
  # we sort by p value follow by ratio however for plot we take top ten and resort by generatio to make it look nice. 
  dfg = rbind ( head ( dfg_full[ dfg_full$group == "up", c("Pathway","p","GeneRatio", "group")], thres ) , head ( dfg_full[ dfg_full$group == "down" , c("Pathway","p","GeneRatio", "group")], thres) )
  
  
  dfg = dfg[ order (dfg$group,   abs ( dfg$GeneRatio) , dfg$p, decreasing = F), ]
  dfg$pathway = factor ( dfg$Pathway, levels =  dfg$Pathway)
  dfg$group = factor ( dfg$group , levels = c("up","down"))  
  
  g1 = ggplot(data=dfg, aes(x=GeneRatio, y=pathway, fill=group)) +
    geom_bar(stat="identity") + scale_fill_manual(values=c(  up="#8bb4c9", down="#f2c830" ), name="NES direction") +
    theme(legend.position="right",  legend.key = element_blank(),
          
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(angle = 90, size=12),
          axis.title.x = element_text(size=12),
          
          axis.title.y     = element_text(size=12), 
          legend.text      =element_text(size=12)
    ) + xlab("Gene ratio") + ylab("") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + ggtitle ( paste( "Hypergeometric TOP positive/negative: ", title ) )
  
  
  
  
  # there might be repeats due to geens being both up and down significant 
  # next we want to exctract the actual genes that were enriched 
  dfg_full2 = dfg_full
  
  
  if ( !is.na(goana_db)){
    
    ingene = data.frame ( stringsAsFactors = F )
    
    # loops through each pathways and pulls out all the gene in that pathways and then intersects with DEG. 
    
    for ( p in dfg_full$Pathway) {
      
      
      ga = unique ( goana_db[goana_db$PathwayID == p, ]$GeneID )
      
      inset = intersect ( unique ( c(up) ), ga ) 
      geneID_up = inset 
      geneName_up = ( refdb[ refdb$entrezgene_id %in% inset,  "external_gene_name" ]  )
      geneName_up = geneName_up[ !duplicated ( geneName_up )]
      
      
      inset = intersect ( unique ( c(down) ), ga ) 
      geneID_down = inset 
      geneName_down = ( refdb[ refdb$entrezgene_id %in% inset,  "external_gene_name" ]  )
      geneName_down = geneName_down[ !duplicated ( geneName_down )]
      
      
      
      
      
      ingene = rbind ( ingene, data.frame ( Pathway = p
                                            , Up_regulated= toString(  geneName_up )
                                            , Down_regulated= toString(geneName_down ), stringsAsFactors = F
      ), stringsAsFactors = F 
      )
      
    }
    
    dfg_full2 = merge ( dfg_full, ingene, by="Pathway")
    #row.names ( dfgFULL ) = dfgFULL$Pathway [ duplicated ( dfgFULL$Pathway ) ]
    # subset 
    
    dfg_full2 = dfg_full2[ order ( dfg_full2$p, -dfg_full2$GeneRatio ), ]
    
  }
  
  
  
  
  return ( list ( plot=g1, df=dfg_full2 ) )
  
}




study.enrich2 = function ( enrich_result, limit=100 , topbubble = 10 ){
  # max_show no longer needed because will add everything 
  ### @@ Start GO   
  # prepare pathways for plotting 
  
  k.go = rbind (
    enrich_result$BP
    , enrich_result$MF
    , enrich_result$KEGG
    , enrich_result$HALLMARK
    , enrich_result$REACT
    , enrich_result$BIOC
    
  )
  
  # we need to study this a bit. 
  # pathways with just 1 gene should be removed
  
  k.go = k.go[ k.go$GeneRatio_Up > .02 | k.go$GeneRatio > .02, ] 
  
  colnames ( k.go )[1] = "Term_Description"
  colnames(k.go)[7] = c("lowest_p")
  k.go$Up_regulated = as.character(k.go$Up_regulated  )
  k.go$Down_regulated = as.character( k.go$Down_regulated)
  k.go = k.go[ , c("Term_Description", "Up_regulated", "Down_regulated", "lowest_p", "GeneRatio", "P.Up" , "P.Down"  ) ]
  k.go$ID = k.go$Term_Description
  
  k.go = k.go[ order ( k.go$lowest_p, -k.go$GeneRatio ), ]
  
  top_go = rbind ( 
    
    head ( k.go[ grepl ( "^BP", k.go$Term_Description) ,], topbubble )
    ,head ( k.go[ grepl ( "^MF", k.go$Term_Description) ,], topbubble )
    ,head ( k.go[ grepl ( "^KEGG", k.go$Term_Description) ,], topbubble )
    ,head ( k.go[ grepl ( "^HALLMARK", k.go$Term_Description) ,], topbubble )
    ,head ( k.go[ grepl ( "^REAC", k.go$Term_Description) ,], topbubble )
  )  
  
  k.go = head ( k.go, limit )
  # organized the terms 
  # 
  k.go.org = cluster_enriched_terms(k.go, plot_clusters_graph = F, method = "hierarchical", use_names=T )
  
  
  
  k.go.repr = k.go.org[k.go.org$Status == "Representative", ]$ID
  
  
  buildkapp.go = data.frame ( stringsAsFactors = F )
  for ( repr in k.go.repr ) {
    
    mem = k.go.org[ k.go.org$ID == repr, ]$Cluster
    mem = k.go.org[k.go.org$Cluster == mem, ]$ID
    mem = mem[ mem != repr ]
    p = k.go.org[ k.go.org$ID == repr, ]$lowest_p
    buildkapp.go = rbind ( buildkapp.go, data.frame ( main = repr, members = toString( paste ( mem , collapse = " , ") ), total = length ( mem), p.value=p  
                                                      , stringsAsFactors = F ) )
    
    
  }
  
  buildkapp.go = buildkapp.go[order (-buildkapp.go$total, buildkapp.go$p ), ]
  buildkapp.go = buildkapp.go 
  
  
  
  # produce plots 
  
  
  
  temp.go =  k.go[ k.go$Term_Description %in% head ( buildkapp.go[buildkapp.go$total >  0, ]$main, 50) , ] 
  temp.go$ID = stringr::str_trunc(temp.go$ID , 45) 
  
  # legacy stuff adopated from gsea 
  # enrich_result$all$all_df = rbind (
  #  enrich_result$BP
  #  , enrich_result$MF
  #  , enrich_result$KEGG
  # , enrich_result$HALLMARK
  #, enrich_result$REACT
  #, enrich_result$BIOC
  
  #) 
  
  plots.go = plot.enrich.post.hyper ( enrich_result = enrich_result , group.go =  buildkapp.go, title="GO", folddf=folddf, maxpath=20)
  
  
  
  
  
  
  # plot first one
  
  
  
  
  
  
  ########### @ test code to plot single network. Please ignore 
  
  if ( 1 > 2 ){
    
    library(igraph)
    ig = c()
    repr = buildkapp.path[1:1, ]$main 
    mem = unlist(strsplit( buildkapp.path[1:1, ]$members, ",")) 
    mem = gsub ( " ","", mem )
    
    
    for ( m in mem ){
      ig = c(ig, c(repr, m ))
      
    }
    
    node.size<-setNames( rep ( 2, length ( mem )) ,c(mem))
    
    
    ig1 <- graph(edges=ig, directed=T)
    lay = layout.reingold.tilford(ig1) 
    
    
    plot(ig1, layout=-lay[, 2:1], vertex.size=node.size,  vertex.label.cex = .8 )
    
    
    
    
    
  }
  
  ####################
  
  # another test plot that is not ncessary. These are bubble plots 
  
  
  ##
  
  return = list ( 
    buildkapp.go =buildkapp.go
    ,plots.go = plots.go
    ,top_go=top_go
  )
  
  
  
  
}

plot.enrich.post <- function ( enrich_result, stat.go , group.go, title="GO") {
  
  geneplot.go = term_gene_graph (  stat.go ) + 
    ggtitle ( paste ( "TOP", title, "Terms"))
  upset.go=      pathfindR::UpSet_plot ( stat.go )   + 
    ggtitle (  paste ( "TOP", title, "Terms"))
  
  
  upset.go.bar =      pathfindR::UpSet_plot ( stat.go , method="barplot")   + 
    ggtitle (  paste ( "TOP", title, "Terms"))
  
  
  
  upset.go = upset.go + geom_vline(xintercept= seq_along(unique(upset.go$data$Term)) + 0.5, color = "grey") 
  upset.go$theme$panel.grid.major.y = element_line(colour="grey", size=0.1)
  
  # change color schema 
  col = c("#99C7E0","#F99356") 
  names ( col ) = c("down","up")
  
  upset.go = upset.go + scale_fill_manual(values = col [ levels ( upset.go$data$Up_Down  ) ]  )
  #upset.path = upset.path + scale_fill_manual(values = col [ levels ( upset.path$data$Up_Down  ) ] )
  
  #geneplot.go + scale_colour_manual(values = c("#989898","#F99356", "#99C7E0" ), name = "", labels = c("term", "DOWN", "UP"))  
  
  # get ggplot bubble 
  
  bubble.go = enrich_result$all[enrich_result$all$Pathway %in% group.go$main[1:10],  ]
  bubble.go = merge ( bubble.go,group.go, by.x="Pathway", by.y="main" )
  bubble.go = bubble.go[ order ( bubble.go$total), ]
  bubble.go$Pathway = stringr::str_trunc(bubble.go$Pathway , 45) 
  bubble.go$Pathway = factor ( bubble.go$Pathway, levels=bubble.go$Pathway)
  #bubble.go$total = factor ( bubble.go$total, levels=bubble.go$total)
  bubble.go$lp = -log2 ( bubble.go$p)
  
  
  bubble.go = ggplot(bubble.go, aes(x=total, y=Pathway ) ) +
    geom_point(aes(color = lp ), size=10) +
    scale_color_gradient2(low="#b5d1eb", high="#ed6218", mid = "#1886ed") + scale_size_area(max_size=5) +
    theme(legend.position="right",  legend.key = element_blank(),
          
          axis.text.y = element_text(size=10),
          axis.text.x = element_text(angle = 90, size=11.5),
          axis.title.x = element_text(size=10),
          
          axis.title.y     = element_text(size=10), 
          legend.text      =element_text(size=10)
    ) + xlab("Total pathway members") + ylab("Main representative pathways") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(color="-log(fdr)") + ggtitle (  paste ( "TOP", title, "Terms") )
  
  return ( list( bubble= bubble.go, upset=upset.go, geneplot=geneplot.go, upset.go.bar=upset.go.bar  ))
}





plot.enrich.post.hyper <- function ( enrich_result , group.go, title="Top Pathways", folddf , maxpath = 20 ) {
  
  
  geneplot.go="none"
  
  temp =  rbind (
    enrich_result$BP
    , enrich_result$MF
    , enrich_result$KEGG
    , enrich_result$HALLMARK
    , enrich_result$REACT
    , enrich_result$BIOC
    
  ) 
  
  bubble.go = temp [temp$Pathway  %in% head ( group.go$main , maxpath ),  ]
  bubble.go$pathway = bubble.go$Pathway
  bubble.go = merge ( bubble.go,group.go, by.x="pathway", by.y="main" )
  bubble.go = bubble.go[ order ( bubble.go$total), ]
  bubble.go$Pathway = gsub ( "_", " ", bubble.go$pathway )
  bubble.go$Pathway = stringr::str_trunc(bubble.go$Pathway , 88) 
  bubble.go$Pathway = factor ( bubble.go$Pathway, levels=bubble.go$Pathway)
  #bubble.go$total = factor ( bubble.go$total, levels=bubble.go$total)
  
  bubble.go$lp = -log2 ( bubble.go$p )
  
  
  bubble.go = ggplot(bubble.go, aes(x=total, y=Pathway ) ) +
    geom_point(aes(color = lp ), size=10) +
    scale_color_gradient2(low="#b5d1eb", high="#ed6218", mid = "#1886ed") + scale_size_area(max_size=5) +
    theme(legend.position="right",  legend.key = element_blank(),
          
          axis.text.y = element_text(size=11),
          axis.text.x = element_text(angle = 90, size=15),
          axis.title.x = element_text(size=20),
          
          axis.title.y     = element_text(size=20), 
          legend.text      =element_text(size=10)
    ) + xlab("Total pathway members") + ylab("Main representative pathways") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(color="-log(fdr)") + ggtitle (  paste ( "TOP", title, "Terms") ) +
    scale_y_discrete(labels = scales:: wrap_format(50))
  
  return ( bubble.go   )
}

plot.enrich.post.cont2 <- function ( enrich_result, stat.go , group.go, title="GO", folddf , maxpath = 20 ) {
  
 
  
  #switch 1 pathway 
  stat.go[1:1, ]$Down_regulated = stat.go[1:1, ]$Up_regulated    
  stat.go[1:1, ]$Up_regulated = ""
  
 
  geneplot.go="none"
  
  upset.go=      pathfindR::UpSet_plot ( stat.go, genes_df =   folddf )   + 
    ggtitle (  paste ( "TOP", title, "Terms"))
  
  
  
  upset.go = upset.go + geom_vline(xintercept= seq_along(unique(upset.go$data$Term)) + 0.5, color = "grey") 
  upset.go$theme$panel.grid.major.y = element_line(colour="grey", size=0.1)
  
  # change color schema 
  #col = c("#99C7E0","#F99356") 
  #names ( col ) = c("down","up")
  
  upset.go = upset.go + scale_fill_gradient( low = "#99C7E0", high = "#F99356")  
  #upset.path = upset.path + scale_fill_manual(values = col [ levels ( upset.path$data$Up_Down  ) ] )
  
  # fix na 
  # there is a strange error where some genes just would not show up and would be plotted as NA so this is a small hack
  na = as.character ( upset.go$data[ is.na ( upset.go$data$Value), ]$Gene ) 
  for ( n in na ){
    value = folddf[folddf$Gene.symbol == n, ]$logFC 
    upset.go$data[ upset.go$data$Gene == n, ]$Value = value 
  }
  
  
  
  #geneplot.go + scale_colour_manual(values = c("#989898","#F99356", "#99C7E0" ), name = "", labels = c("term", "DOWN", "UP"))  
  
  # get ggplot bubble 
  
  bubble.go = enrich_result$all[enrich_result$all_df$pathway %in% head ( group.go$main , maxpath ),  ]
  bubble.go = merge ( bubble.go,group.go, by.x="pathway", by.y="main" )
  bubble.go = bubble.go[ order ( bubble.go$total), ]
  bubble.go$Pathway = gsub ( "_", " ", bubble.go$pathway )
  bubble.go$Pathway = stringr::str_trunc(bubble.go$Pathway , 88) 
  bubble.go$Pathway = factor ( bubble.go$Pathway, levels=bubble.go$Pathway)
  #bubble.go$total = factor ( bubble.go$total, levels=bubble.go$total)
  
  bubble.go$lp = -log2 ( bubble.go$padj)
  
  
  bubble.go = ggplot(bubble.go, aes(x=total, y=Pathway ) ) +
    geom_point(aes(color = lp ), size=10) +
    scale_color_gradient2(low="#b5d1eb", high="#ed6218", mid = "#1886ed") + scale_size_area(max_size=5) +
    theme(legend.position="right",  legend.key = element_blank(),
          
          axis.text.y = element_text(size=11),
          axis.text.x = element_text(angle = 90, size=15),
          axis.title.x = element_text(size=20),
          
          axis.title.y     = element_text(size=20), 
          legend.text      =element_text(size=10)
    ) + xlab("Total pathway members") + ylab("Main representative pathways") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(color="-log(fdr)") + ggtitle (  paste ( "TOP", title, "Terms") ) +
    scale_y_discrete(labels = scales:: wrap_format(50))
  
  return ( list( bubble= bubble.go, upset=upset.go, geneplot=geneplot.go ))
}

plot.enrich.post.cont <- function ( enrich_result, stat.go , group.go, title="GO", folddf ) {
  
  
  
  #switch 1 pathway 
  stat.go[1:1, ]$Down_regulated = stat.go[1:1, ]$Up_regulated
  stat.go[1:1, ]$Up_regulated = ""
  
  
  
 
  geneplot.go="none"
  
  upset.go=      pathfindR::UpSet_plot ( stat.go, genes_df =   folddf )   + 
    ggtitle (  paste ( "TOP", title, "Terms"))
  
  
  
  upset.go = upset.go + geom_vline(xintercept= seq_along(unique(upset.go$data$Term)) + 0.5, color = "grey") 
  upset.go$theme$panel.grid.major.y = element_line(colour="grey", size=0.1)
  
  # change color schema 
  #col = c("#99C7E0","#F99356") 
  #names ( col ) = c("down","up")
  
  upset.go = upset.go + scale_fill_gradient( low = "#99C7E0", high = "#F99356")  
  #upset.path = upset.path + scale_fill_manual(values = col [ levels ( upset.path$data$Up_Down  ) ] )
  
  # fix na 
  # there is a strange error where some genes just would not show up and would be plotted as NA so this is a small hack
  na = as.character ( upset.go$data[ is.na ( upset.go$data$Value), ]$Gene ) 
  for ( n in na ){
    value = folddf[folddf$Gene.symbol == n, ]$logFC 
    upset.go$data[ upset.go$data$Gene == n, ]$Value = value 
  }
  
  
  
  #geneplot.go + scale_colour_manual(values = c("#989898","#F99356", "#99C7E0" ), name = "", labels = c("term", "DOWN", "UP"))  
  
  # get ggplot bubble 
  
  bubble.go = enrich_result$all[enrich_result$all$Pathway %in% group.go$main[1:10],  ]
  bubble.go = merge ( bubble.go,group.go, by.x="Pathway", by.y="main" )
  bubble.go = bubble.go[ order ( bubble.go$total), ]
  bubble.go$Pathway = stringr::str_trunc(bubble.go$Pathway , 45) 
  bubble.go$Pathway = factor ( bubble.go$Pathway, levels=bubble.go$Pathway)
  #bubble.go$total = factor ( bubble.go$total, levels=bubble.go$total)
  bubble.go$lp = -log2 ( bubble.go$pval)
  
  
  bubble.go = ggplot(bubble.go, aes(x=total, y=Pathway ) ) +
    geom_point(aes(color = lp ), size=10) +
    scale_color_gradient2(low="#b5d1eb", high="#ed6218", mid = "#1886ed") + scale_size_area(max_size=5) +
    theme(legend.position="right",  legend.key = element_blank(),
          
          axis.text.y = element_text(size=10),
          axis.text.x = element_text(angle = 90, size=11.5),
          axis.title.x = element_text(size=10),
          
          axis.title.y     = element_text(size=10), 
          legend.text      =element_text(size=10)
    ) + xlab("Total pathway members") + ylab("Main representative pathways") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(color="-log(fdr)") + ggtitle (  paste ( "TOP", title, "Terms") )
  
  return ( list( bubble= bubble.go, upset=upset.go, geneplot=geneplot.go ))
}




exp_cells <- function ( gene_query="", sob , assay_this = "SCT", mthis="seurat_clusters"
                        , min=20
                        , qprob=.25, 
                        avethis = "all" # set to all or nonzero to take average of all or just non-zero
                        ,minsize= 2.5, maxsize=8, rorder =F, scaleby = "zscore", xlab_size = 15, ylab_size = 15
){
  
  cells = sob[[assay_this]]$data 
  gene_indices <- which(rownames(cells) %in% gene_query )
  
  cells = reshape2::melt ( as.matrix  ( cells[gene_indices,  , drop=F]) )
  colnames ( cells ) = c("gene","cell", "value")
  
  cells = merge ( cells, sob@meta.data[ , mthis, drop=F], by.x="cell", by.y="row.names")
  colnames(cells)[ grepl (mthis, colnames(cells))] = "bin"
  dim ( cells )
  # calculate the average non-zero expression and how percent hits
  q = as.numeric ( quantile ( cells[cells$value != 0, ] $value, qprob) )
  
  result <- cells %>%
    group_by(gene,bin) %>%
    # Compute total cells and non-zero cells for each gene
    summarise(
      total_cells = n_distinct(cell),
      non_zero_cells = n_distinct(cell[value != 0]),
      q25_cells = n_distinct(cell[value > q]),
      avg_value = mean( log1p ( value[value != 0]) ),
      median_value = median(value[value != 0]),
      min_value = min(value[value != 0]),
      max_value = max(value[value != 0]), 
      ave_value_all =  mean ( log1p (value) ) , 
      q25 = as.numeric ( quantile( value[value != 0], prob = qprob ) )
    ) %>%
    # Compute the percentage of cells with non-zero values
    mutate(percent_cells_non_zero = (non_zero_cells / total_cells) * 100
           ,percent_q25_cells = (q25_cells / total_cells) * 100
    ) %>%
    select(bin, gene, percent_cells_non_zero, avg_value, median_value, min_value, max_value, ave_value_all, percent_q25_cells, q25 )
  
  result <-  data.frame ( result , stringsAsFactors = F)
  
  # Replace -Inf, Inf, and NaN with 0 for all numeric columns
  numeric_cols <- sapply(result, is.numeric)
  result[numeric_cols] <- lapply(result[numeric_cols], 
                                 function(x) replace(x, is.infinite(x) | is.nan(x) | is.na(x), 0))
  
  result = result[ order ( result$bin), ]
  result$bin = factor ( result$bin, levels = sort ( unique ( result$bin)))
  result$label = result$percent_cells_non_zero
  result$label25 = result$percent_q25_cells
  
  
  result$label[result$label < min] <- 0
  result$label25[result$label < min] <- 0
  result$gene = factor ( result$gene, levels = rev ( gene_query [ gene_query %in% unique ( cells$gene)]
  ))
  
  if ( avethis == "all"){
    result$avethis = result$ave_value_all
  } else {
    result$avethis = result$avg_value
  }
  
  if ( rorder == TRUE ){
    result = result [ order ( result$avethis, decreasing = T), ]
    result$bin = factor ( result$bin, levels = unique ( as.character ( result$bin) ))
  }
  
  # scale 
  minmax <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }
  
  z_score <- function(x) {
    (x - mean(x)) / sd(x)
  }
  
  
  if ( scaleby == "zscore"){    
    result <- result %>%
      group_by(gene) %>%
      mutate(scaled_avethis = z_score(avethis))
  } else if ( scaleby=="minmax" )  {
    result <- result %>%
      group_by(gene) %>%
      mutate(scaled_avethis = minmax(avethis))
  } else if ( scaleby=="log2" )  {
    result <- result %>%
      group_by(gene) %>%
      mutate(scaled_avethis = log2(avethis+1))
  } else {
    result$scaled_avethis = result$avethis
  }
  
  
  
  plot_25 =  ggplot(result, aes(x = bin, y = gene, size = ifelse(label25 == 0, NA, label25), color = scaled_avethis )) +
    geom_point(alpha = 0.6) +
    #scale_color_gradient2(low = "lightgrey",   high = "#b8241a")  +
    scale_colour_gradientn(colours = c("lightgrey","grey", "orange", "#fca503","#ad6661", "#b8241a")) +
    labs(title = "Bubble Plot",
         x = mthis,
         y = "",
         size = "% Cells Non-Zero",
         color = "Average Exp") +
    theme_minimal() + scale_size(range = c(minsize, maxsize)) +
    theme(panel.grid.major = element_blank()
          , panel.grid.minor = element_blank()
          ,panel.background = element_blank()
          , axis.line = element_line(colour = "black")
    ) +
    theme(legend.position="right",  legend.key = element_blank(),
          # element_blank()
          axis.text.y = element_text(size= ylab_size ),
          axis.text.x =  element_text(size= xlab_size , angle = 45, hjust = 1 ),
          axis.title.x = element_text(size=25),
          axis.title.y     = element_text(size=25), 
          legend.text      =element_text(size=25),
          legend.title = element_text(size=25)
          #legend.key.size = unit(2, 'cm'),
          #plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
          # hjust centers the title
    ) + ggtitle ( "")+ xlab ("Cluster Assignment")
  
  plot_0 =  ggplot(result, aes(x = bin, y = gene, 
                               size = ifelse(label == 0, NA, label)
                               , color =  scaled_avethis  )) +
    geom_point(alpha = 0.6) +
    #scale_color_gradient2(low = "lightgrey",   high = "#b8241a")  +
    scale_colour_gradientn(colours = c("lightgrey","grey", "orange", "#fca503", "#ad6661", "#b8241a")) +
    labs(title = "Bubble Plot",
         x = mthis,
         y = "",
         size = "% Cells Non-Zero",
         color = "Average Exp") +
    theme_minimal() + scale_size(range = c(minsize, maxsize)) +
    theme(panel.grid.major = element_blank()
          , panel.grid.minor = element_blank()
          ,panel.background = element_blank()
          , axis.line = element_line(colour = "black")
    ) +
    theme(legend.position="right",  legend.key = element_blank(),
          # element_blank()
          axis.text.y = element_text(size= ylab_size ),
          axis.text.x =  element_text(size= xlab_size , angle = 45, hjust = 1 ),
          axis.title.x = element_text(size=25),
          axis.title.y     = element_text(size=25), 
          legend.text      =element_text(size=25),
          legend.title = element_text(size=25)
          #legend.key.size = unit(2, 'cm'),
          #plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
          # hjust centers the title
    ) + ggtitle ( "") + xlab ("Cluster Assignment")
  
  return ( list ( result=result, plot_0 = plot_0, plot_25 = plot_25))
  
}


# summarize qc 
summqc <- function ( qcd, feature, n){
  summary_vector <- unlist(summary(qcd))
  summary_df <- data.frame(matrix(summary_vector, ncol = length(summary_vector), byrow = TRUE))
  colnames(summary_df) <- names(summary_vector)
  summary_df$sample = n 
  summary_df$feature = feature
  summary_df = summary_df[ , c("sample",  "feature", names(summary_vector) )]
  return ( summary_df )
}



plot_ratios <- function ( sob, cmp_this, x_this, ccc_temp , sqrt_yes = 0  ){
  
  sob[["cmp_this"]] = sob[[cmp_this]]
  sob[["x_this"]] = sob[[x_this]]
  
  sobmeta = sob@meta.data
  
  total_pid = data.frame ( table ( sobmeta$orig.ident ) )
  total_pid = setNames(total_pid$Freq, as.character ( total_pid$Var1))
  
  
  total_cellpid = sobmeta %>% group_by( x_this, orig.ident, cmp_this ) %>% 
    summarise(Total=n()) %>% data.frame( )
  
  colnames ( total_cellpid ) [which ( colnames ( total_cellpid ) == "orig.ident" )] = "pid"
  # multiply by 100 for better visual 
  
  total_cellpid$percent = total_cellpid$Total / total_pid[ total_cellpid$pid] *100
  
  
  # sanity check   
  total_cellpid %>% group_by( pid ) %>% summarise(total=sum(percent))
  total_cellpid %>% group_by( x_this, cmp_this ) %>% summarise(mean=mean(percent)) %>% data.frame()
  
  
  # remove unknown clusters and celltype
  total_cellpid=total_cellpid[! ( total_cellpid[, "x_this"] %in% c ( "Unclassified", "Unknown", "undetermined") |
                                    total_cellpid[, "cmp_this"] %in% c ( "Unclassified", "Unknown", "undetermined")  
                                  
  )
  , ]
  
  
  ylab = "cell %"
  if ( sqrt_yes == 1){
    total_cellpid$sqrt = sqrt(total_cellpid$percent )
    ylab = "sqrt ( cell % ) "
  }else {
    total_cellpid$sqrt = c(total_cellpid$percent )
  }
  
  
  count_specp <- total_cellpid %>%
    group_by(x_this, cmp_this) %>%
    summarise(mean_norm = mean(sqrt), se_norm = sd(sqrt)/sqrt(n()) ) %>% data.frame()
  
  count_specp$x_this = as.character(count_specp$x_this)
  
  ranked_data <- count_specp %>%
    arrange(desc(cmp_this), desc(mean_norm))
  
  count_specp$x_this = factor(as.character ( count_specp$x_this), level=unique ( ranked_data$x_this) )
  
  
  pv_result <- total_cellpid  %>%
    group_by(x_this) %>%
    nest() %>% # important to treat each celltype as a dataframe
    mutate(test_result = map(data, ~wilcox.test( sqrt ( percent ) ~ cmp_this, data = .)),
           p_value = map_dbl(test_result, ~.$p.value),
           statistic = map_dbl(test_result, ~.$statistic) ) %>%
    ungroup() %>%
    mutate(FDR = p.adjust(p_value, method = "fdr")
           ,bonferroni = p.adjust(p_value, method = "bonferroni")
    )
  
  pv_result$data = NULL 
  pv_result$test_result = NULL 
  pv_result$statistic = NULL 
  pv_result$bonferroni = NULL # only use this if you have extremely low p value 
  
  
  pv_result = pv_result %>%
    mutate(pv = case_when(
      p_value < 0.01 & FDR < 0.05 ~ "**",
      p_value < 0.05 & p_value >= 0.01 & FDR < 0.05 ~ "*",
      TRUE ~ ""
    ))
  
  pv_result = pv_result %>% mutate_if(is.numeric, round, digits = 3 ) %>% data.frame()
  
  # calculate fold change 
  
  
  fold_change_df <- total_cellpid %>%
    group_by(x_this, cmp_this) %>%
    summarise(mean_percent = mean(percent)) %>%
    spread(key = cmp_this, value = mean_percent) %>%
    mutate(log2FC = round ( log2 ( jmml / healthy_cntr ), 2)
    )
  
  pv_result = merge (fold_change_df, pv_result, by= "x_this" )
  
  
  # Compute the max mean_norm for each CellType
  max_values <- count_specp %>% 
    group_by(x_this) %>%
    summarise(max_mean_norm = max(mean_norm + se_norm))
  
  # Merge the computed max values and pv into count_specp
  count_specp2 <- left_join(count_specp, max_values, by = "x_this")
  count_specp2 <- left_join(count_specp2, pv_result[, c("x_this", "pv")], by = "x_this")
  count_specp2$x_this = factor ( count_specp2$x_this, level = unique ( ranked_data$x_this))
  
  
  
  bar_celltype = ggplot(count_specp2, aes(x = x_this, y = mean_norm, fill = cmp_this)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = mean_norm - se_norm, ymax = mean_norm + se_norm), 
                  width = 0.2, position = position_dodge(0.9), color = "#8f8f8f" ) +
    
    
    #geom_text(data = distinct(count_specp2, CellType, .keep_all = TRUE), size = 7,
    #         aes(label = ifelse(pv == "**", " ** ", 
    #                           ifelse(pv == "*", " * ", "")), 
    #           y = max_mean_norm + 0.51), vjust = 0) +
    
    ylab("Norm") +
    xlab("Celltype") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ylab(ylab) + 
    xlab("") + 
    theme(legend.position="bottom",  legend.key = element_blank(),
          axis.text.y = element_text(size= 15 ),
          axis.text.x = element_text(angle = 45, size= 11 , hjust=1 ),
          axis.title.x = element_text(size=12),
          axis.title.y     = element_text(size=15), 
          legend.text      =element_text(size=12),
          plot.title = element_text(size = 15, face = "bold")
    )  +  scale_fill_manual(values= ccc_temp  ) 
  
  
  bar_celltype
  
  
  bar_celltype_noerror = ggplot(count_specp2, aes(x = x_this, y = mean_norm, fill = cmp_this)) +
    geom_bar(stat = "identity", position = "dodge") +
    
    ylab("Norm") +
    xlab("Celltype") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    ylab(ylab) + 
    xlab("") + 
    theme(legend.position="bottom",  legend.key = element_blank(),
          axis.text.y = element_text(size= 15 ),
          axis.text.x = element_text(angle = 45, size= 11 , hjust=1 ),
          axis.title.x = element_text(size=12),
          axis.title.y     = element_text(size=15), 
          legend.text      =element_text(size=12),
          plot.title = element_text(size = 15, face = "bold")
    )  +  scale_fill_manual(values= ccc_temp  ) 
  
  
  
  pv_result = pv_result[match(unique ( ranked_data$x_this) , pv_result$x_this),]
  pv_result = pv_result[!is.na ( pv_result$p_value), ]
  
  pretty_stat <- ggtexttable(pv_result, rows = NULL,  theme = ttheme("light"
                                                                     #, base_size = 18 
  )) 
  pretty_stat <- table_cell_font (pretty_stat, column = 4, row= c(1:nrow ( pv_result)+1),
                                  face = "bold", size=15) %>% 
    tab_add_vline(at.column = 2, column.side = "left", from.row = 2, linetype = 2) %>% 
    tab_add_vline(at.column = 4, column.side = "left", from.row = 2, linetype = 2) 
  
  
  ## put it all together 
  library(patchwork)
  layout <- "
AAAAACCCCC
AAAAACCCCC
AAAAADDDDD
BBBBBDDDDD
BBBBBDDDDD
BBBBBDDDDD
BBBBBDDDDD
"
  pretty_stat2=pretty_stat
  
  #cell_type_plot_man + sting1 + bar_celltype + pretty_stat2  + 
  # plot_layout(design = layout)   + plot_annotation(tag_levels = 'A')
  
  return ( list (
    bar_celltype=bar_celltype,
    pretty_stat=pretty_stat,
    bar_celltype_noerror = bar_celltype_noerror
    
  ))
  
}





## set colors

# generate colors 
getPalette = colorRampPalette(brewer.pal(9, "Set1")) # expand color pallete
main_clrs = readRDS( config$precolor )

rcolor = main_clrs$rcolor
main_clrs = main_clrs$main_clrs


group_color = c(  jmml = "#B89977", healthy_cntr = "#97BD5E")

ctccc <- c("HSC", "mDCs", "C_Monocytes", "CD8+ NKT", "GMP", "Pre-B", "N_Monocytes", "Erythroid", "Neutrophils", 
           "Naive CD4+ T", "MLP", "Pro-B", "ISG+", "Naive B", "MΦ")


ccc_celltype = setNames( main_clrs[ 1:length ( unique ( ctccc ) )]  
                         , unique ( ctccc)
)


ccc_celltype [ "Naive CD4+ T"] = "grey"
ccc_celltype["Erythroid"]= "#c75b54"
ccc_celltype["ISG+"]= "red"
ccc_celltype["mDCs"] = "#9964e3"
ccc_celltype["C_Monocytes" ] = "#96795d"
ccc_celltype["GMP"] = "#539e49"
ccc_celltype["MLP"] = "#1ff2cc"
ccc_celltype[ "Neutrophils"] = "#edbc0c"
ccc_celltype[ "MΦ" ] = "#ab18a6"


# simple themes for multiple UMAP plots 

umapSimple_theme <- function() {

  list(
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      #axis.line = element_blank(),
      plot.title = element_text(size = 15, face = "bold")
    ),
    labs( x = "", y = "", title = "")
  )
  
  
}

## cleanup split maps 

cleanupsplits <- function ( xplot){
  xplot[[1]] = xplot[[1]] + umapSimple_theme() + ggtitle ( "Healthy Controls")
  xplot[[2]] = xplot[[2]] + umapSimple_theme() + ggtitle ( "JMML")
  return ( xplot )
}


