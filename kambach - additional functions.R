###########################################
# Additional Functions for FeedBaCks
# script by: Dr. Stephan Kambach (stephan.kambach@gmail.com)
# initial date: 18.03.2023
###########################################
fetch.european.borders = function(){
  
  some.eu.countries <- c(
    "Algeria", "Austria", "Azerbajan", "Belgium", "Bulgaria", "Croatia",
    "Cyprus", "Czech republic", "Denmark", "Egypt", "Estonia",
    "Finland", "France", "Georgia",  "Germany", "Greece", 
    "Greenland", "Hungary", "Iceland", "Iran", "Ireland", "Italy", "Kazakhstan", 
    "Latvia", "Libya", "Lithuania", "Luxembourg", "Malta", "Morocco", 
    "Netherlands", "Norway", "Poland", "Portugal", "Romania", "Russia", "Slovakia",
    "Slovenia", "Spain", "Sweden", "Syria", "Tunisia", "Turkey", "UK", "Ukraine")
  
  # Retrievethe map data
  some.eu.maps <- map_data("world") %>%
    rename(longitude = long, latitude = lat)
  
  # Compute the centroid as the mean longitude and lattitude
  # Used as label coordinate for country's names
  region.lab.data <- some.eu.maps %>%
    group_by(region) %>%
    summarise(longitude = mean(longitude), latitude = mean(latitude))
  
  results = list("borders" = some.eu.maps,
                 "labels" = region.lab.data)
  return(results)
}

#################################
# read.in.EVA.header.data = function(){
#   results = read_delim("data/eva/current working version/FeedBacks20210512_corrected_header.csv", delim = "\t",
#                              col_types = "iiccDdiidddddcicccccdiccci") %>%
#     rename(PlotObservationID = "PlotObservationID",
#            plotID = "PlotID",
#            plotnr = "TV2 relev? number",
#            country = "Country",
#            ref = "Biblioreference",
#            ref_tab = "Nr. table in publ.",
#            ref_tab_nr = "Nr. relev? in table",
#            scale = "Cover abundance scale",
#            project = "Project",
#            author = "Author",
#            date = "Date of recording",
#            syntaxon = "Syntaxon",
#            plot_area = "Relev? area (m?)",
#            utm_code = "UTM grid system code",
#            plot_elevation = "Altitude",
#            plot_aspect = "Aspect (?)",
#            plot_slope = "Slope (?)",
#            cover_total = "Cover total (%)",
#            cover_trees = "Cover tree layer (%)",
#            cover_shrubs = "Cover shrub layer (%)",
#            cover_herbs = "Cover herb layer (%)",
#            cover_mosses = "Cover moss layer (%)",
#            cover_lichens = "Cover lichen layer (%)",
#            cover_algae = "Cover algae layer (%)",
#            cover_litter = "Cover litter layer (%)",
#            cover_water = "Cover open water (%)",
#            cover_rock = "Cover bare rock (%)",
#            height_tree_largest = "Height (highest) trees (m)",
#            height_tree_lowest = "Height lowest trees (m)",
#            height_shrub_heighest = "Height (highest) shrubs (m)",
#            height_shrub_lowest = "Height lowest shrubs (m)",
#            height_herbs_highest = "Aver. height (high) herbs (cm)",
#            height_herbs_lowest = "Aver. height lowest herbs (cm)",
#            height_herbs_max = "Maximum height herbs (cm)",
#            height_cryptogams_max = "Maximum height cryptogams (mm)",
#            identified_mosses = "Mosses identified (Y/N)",
#            identified_lichens = "Lichens identified (Y/N)",
#            locality = "Locality",
#            name_association = "Name association",              
#            name_alliance = "Name alliance",
#            author_nam = "AUTHOR_NAM",
#            expert_system = "Expert System",
#            longitude = "Longitude",
#            latitude = "Latitude",
#            location_uncertainty = "Location uncertainty (m)",
#            dataset = "Dataset",
#            access = "Access regime") %>%
#     mutate(date = as_date(date, format = "%d.%m.%Y")) %>%
#     mutate(year = year(date))
#   
#   return(results)
# }
# 
###############################
create.raster.tibble = function(dat.all, ncols, nrows, field){
  
  # create dat.temp to restrict plot to actual data
  dat.temp = dat.all %>%
    drop_na(all_of(field))
  
  # points into a SpatialPointsDataFrame object
  coords <- cbind(x=dat.temp[["longitude"]],y=dat.temp[["latitude"]])
  sPDF <- SpatialPointsDataFrame(coords,data=dat.temp)

  # set extend
  xmn <- min(dat.temp[["longitude"]]) 
  xmx <- max(dat.temp[["longitude"]]) 
  ymn <- min(dat.temp[["latitude"]]) 
  ymx <- max(dat.temp[["latitude"]]) 
  
  # create grid
  blankRaster <- raster(nrows=nrows, ncols=ncols, xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx)
  
  #adding data into raster to avoid 'no data' error
  blankRaster[] <- 0
  
  # calculate number of observations per grid cell
  blankRaster[] <- 0
  rastercountsPoints <- rasterize(x=sPDF, y=blankRaster, field= field , fun="count")
  counts = tibble(as.data.frame(rastercountsPoints, xy = T, na.rm = T)) %>%
    rename(longitude = x, latitude = y, value = layer)
  
  # calculate mean (or other function) of points per cell 
  blankRaster[] <- 0
  rasterMeanPoints <- rasterize(x=sPDF, y=blankRaster, field= field , fun=mean)
  values = tibble(as.data.frame(rasterMeanPoints, xy = T, na.rm = T)) %>%
    rename(longitude = x, latitude = y, value = layer)
  
  results = list(values = values,
                 counts = counts)
  
  return(results)
}

######################
# Calculate CWV
variance2.fun <- function(trait, abu){
  res <- as.double(NA)
  abu <- abu[!is.na(trait)]
  trait <- trait[!is.na(trait)]
  abu <- abu/sum(abu)
  if (length(trait)>1){
    # you need more than 1 observation to calculate variance
    # for calculation see 
    # http://r.789695.n4.nabble.com/Weighted-skewness-and-curtosis-td4709956.html
    m.trait <- weighted.mean(trait,abu)
    res <- sum(abu*(trait-m.trait)^2)}
  res}

#######################
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

##########################
plot.my.varpart = function(A, B, C, AB, AC, BC, ABC, Random,
                           A.name, B.name, C.name, title.name,
                           save.file, plot.height, plot.width){
  
  totalexpl <- ABC
  Apure <- totalexpl - BC
  Bpure <- totalexpl - AC
  Cpure <- totalexpl - AB
  ABoverlap <- totalexpl - Apure - Bpure - C
  BCoverlap <- totalexpl - Bpure - Cpure - A
  ACoverlap <- totalexpl - Apure - Cpure - B
  ABCoverlap <- totalexpl - Apure - Bpure - Cpure - ABoverlap - BCoverlap - ACoverlap
  unexpl <- 100 - totalexpl - Random
  
  output.names <- c(paste("Pure", A.name), 
                    paste("Pure", B.name), 
                    paste("Pure", C.name), 
                    paste0("Pure ", A.name, "_", B.name, " overlap"), 
                    paste0("Pure ", B.name, "_", C.name, " overlap"), 
                    paste0("Pure ", A.name, "_", C.name, " overlap"), 
                    paste0(A.name, "_", B.name, "_", C.name, " overlap"),
                    "Random",
                    "Unexplained")
  
  results <- data.frame(c(Apure, Bpure, Cpure, ABoverlap, BCoverlap, ACoverlap, ABCoverlap, Random, unexpl), 
                        row.names = output.names)
  colnames(results) <- "Proportion"
  
  totalexpl <- ifelse(totalexpl >= 0.01, round(totalexpl, 2), "")
  unexpl <- ifelse(unexpl >= 0.01, round(unexpl, 2), "")
  Random = ifelse(Random >= 0.01, round(Random, 2), "")
  Apure <- ifelse(Apure >= 0.01, round(Apure, 2), "")
  Bpure <- ifelse(Bpure >= 0.01, round(Bpure, 2), "")
  ABoverlap <- ifelse(ABoverlap >= 0.01, round(ABoverlap, 2), "")
  Cpure <- ifelse(Cpure >= 0.01, round(Cpure, 2), "")
  BCoverlap <- ifelse(BCoverlap >= 0.01, round(BCoverlap, 2), "")
  ACoverlap <- ifelse(ACoverlap >= 0.01, round(ACoverlap, 2), "")
  ABCoverlap <- ifelse(ABCoverlap >= 0.01, round(ABCoverlap, 2), "")
  
  dat.circles =  data.frame(x0 = c(1,2,3),
                            y0 = c(1,-0.5,1),
                            r = c(1.5,1.5,1.5),
                            fill = c("red", "blue", "green"))
  
  dat.label =  data.frame(
    x = c(0.5, 3.5,  2,    2,  1,  3, 2),
    y = c(1.2, 1.2, -0.8,  1.5,  0,  0, 0.5),
    label = c(Apure, Bpure, Cpure, ABoverlap, ACoverlap, BCoverlap, ABCoverlap)
  )
  
  plot.temp = ggplot() +
    geom_circle(data = dat.circles, 
                aes(x0 = x0, y0 = y0, r = r, fill = fill), alpha = 0.3) +
    geom_rect(aes(xmin = - 1, xmax = 5, ymin = -2.5, ymax = 3), fill = NA, color = "black") + 
    geom_polygon(aes(x = c(-1, -1, 1), y = c(-2.5, 0, -2.5)), alpha = 0.3) + 
    geom_segment(aes(x = -1, y = 0, xend = 1, yend = -2.5)) + 
    geom_text(data = dat.label, 
              aes(x = x, y = y, label = label)) +
    geom_text(aes(x = c(1, 3, 2), y = c(2, 2, -1.5), label = c(A.name, B.name, C.name)), size = 6) +
    geom_text(aes(x = -0.8, y = -2, label = paste0(Random, "\nRandom effect")), color = "black", hjust = 0) +
    geom_text(aes(x = 4.8, y = -2, label = paste0(unexpl, "\nunaccounted")), color = "black", hjust = 1) +
    geom_label(aes(x = -1, y = 3.3, label = title.name), hjust = 0) + 
    theme_void() + 
    theme(panel.grid = element_blank(), legend.position = "none", panel.background = element_rect(fill = "white"))
  
  if (all.equal(sum(results, na.rm = TRUE), 100)){
    cat("")
  }else{
    warning("Results don't sum up to 1; are you sure your input data are correct?")}
  
  ggsave(plot = plot.temp,
         filename = save.file,
         height = plot.height, 
         width = plot.width)
  
  plot(plot.temp)
  return(results)
}


##########################
change.strip.colors = function(plot.temp, colors){
  g <- ggplot_gtable(ggplot_build(plot.temp))
  strip_both <- which(grepl('strip-', g$layout$name))
  fills <- colors
  k <- 1
  for (i in strip_both) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  return(g)
}

###################
get.variable.importance.gam = function(gam.full, dat.full){
  
  results = data.frame("trait" = str_sub(names(gam.full$sp),7,-2),
                       "deviance" = NA)
  
  formula.full = formula(gam.full)
  terms.full = names(gam.full$sp)
  
  # intercept-only model
  formulat.intercept = update(formula.full, paste("~  1"))
  gam.intercept = gam(formulat.intercept, data = dat.full)
  
  for(i in 1:length(terms.full)){
    
    # select term
    term.temp = terms.full[i]
    
    # reduce formula
    formulat.temp = update(formula.full, paste("~ . -", term.temp))
    
    # run nested model
    gam.nested = update(gam.full, formulat.temp,
                        sp = c(gam.full$sp[-i]))
    
    results$deviance[i] = 100 * (deviance(gam.nested)-deviance(gam.full))/deviance(gam.intercept)
    print(paste0("deviance calculated: ", results$trait[i]))
  }
  return(results)
}

#############################
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = length(list(...)), position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

###################################
make.a.nice.pca.plot = function(
    data.temp, nfactors.temp = 8, rotate.temp = "varimax", print.plot = F, recode.trait.names, switch.axes.temp = NULL){
  
  pca.unrotated.temp = psych::principal(r = data.temp, nfactors = nfactors.temp, rotate = "none")
  pca.rotated.temp = psych::principal(r = data.temp, nfactors = nfactors.temp, rotate = rotate.temp)
  
  data.to.plot = list(ind_temp = as_tibble(pca.rotated.temp$scores) %>%
                        setNames(paste("Dim.", 1:nfactors.temp, sep = "")),
                      
                      var_temp = data.frame(unclass(pca.rotated.temp$loadings)) %>%
                        rownames_to_column(var = "trait") %>% 
                        as_tibble() %>%
                        setNames(c("trait", paste("Dim.", 1:nfactors.temp, sep = ""))),
                      var_accounted = data.frame(pca.rotated.temp$Vaccounted) %>%
                        rownames_to_column(var = "metric") %>% 
                        as_tibble() %>%
                        setNames(c("metric", paste("Dim.", 1:nfactors.temp, sep = ""))))
  
  # rotate where necessary
  if(!is_empty(switch.axes.temp)){
    for(i in 1:length(switch.axes.temp)){
      data.to.plot$ind_temp[,switch.axes.temp] = - data.to.plot$ind_temp[,switch.axes.temp]
      data.to.plot$var_temp[,switch.axes.temp] = - data.to.plot$var_temp[,switch.axes.temp]}}
  
  # shorten trait names
  data.to.plot$var_temp$trait = recode.trait.names(data.to.plot$var_temp$trait)
  
  if(nfactors.temp >= 1){
    scale.segment.dim1 = 0.8 * max(abs(quantile(data.to.plot$ind_temp$Dim.1, c(0.001, 0.999)))) / 
    max(abs(data.to.plot$var_temp$Dim.1))}
  if(nfactors.temp >= 2){
    scale.segment.dim2 = 0.8 * max(abs(quantile(data.to.plot$ind_temp$Dim.2, c(0.001, 0.999)))) / 
      max(abs(data.to.plot$var_temp$Dim.2))}
  if(nfactors.temp >= 3){
    scale.segment.dim3 = 0.8 * max(abs(quantile(data.to.plot$ind_temp$Dim.3, c(0.001, 0.999)))) / 
      max(abs(data.to.plot$var_temp$Dim.3))}
  if(nfactors.temp >= 4){
    scale.segment.dim4 = 0.8 * max(abs(quantile(data.to.plot$ind_temp$Dim.4, c(0.001, 0.999)))) / 
      max(abs(data.to.plot$var_temp$Dim.4))}
  if(nfactors.temp >= 5){
    scale.segment.dim5 = 0.8 * max(abs(quantile(data.to.plot$ind_temp$Dim.5, c(0.001, 0.999)))) / 
      max(abs(data.to.plot$var_temp$Dim.5))}
  if(nfactors.temp >= 6){
    scale.segment.dim6 = 0.8 * max(abs(quantile(data.to.plot$ind_temp$Dim.6, c(0.001, 0.999)))) / 
      max(abs(data.to.plot$var_temp$Dim.6))}
  if(nfactors.temp >= 7){
    scale.segment.dim7 = 0.8 * max(abs(quantile(data.to.plot$ind_temp$Dim.7, c(0.001, 0.999)))) / 
      max(abs(data.to.plot$var_temp$Dim.7))}
  if(nfactors.temp >= 8){
    scale.segment.dim8 = 0.8 * max(abs(quantile(data.to.plot$ind_temp$Dim.8, c(0.001, 0.999)))) / 
      max(abs(data.to.plot$var_temp$Dim.8))}
  
  # remove "CWM_"
  data.to.plot$var_temp$trait = gsub("CWM_", "", data.to.plot$var_temp$trait)
  
  plots = list()
  
  if(nfactors.temp >= 2){
    plots[["pca12"]] = ggplot() +
      geom_hex(data = data.to.plot$ind_temp, aes(x = Dim.1, y = Dim.2), 
               bins = 100) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_segment(data = data.to.plot$var_temp, 
                   aes(x = 0, y = 0, xend = Dim.1 * scale.segment.dim1, yend = Dim.2 * scale.segment.dim2),
                   arrow = arrow(length = unit(0.8, "cm"))) +
      geom_label(data = data.to.plot$var_temp, aes(x = Dim.1 * scale.segment.dim1 * 1.2, y = Dim.2 * scale.segment.dim2 * 1.2, label = trait), 
                 alpha = 0.8) +
      scale_fill_gradientn(colours=c("#E8E8E8","#808080"),guide = "none", na.value = NA) +
      xlab(paste0("Axis 1: ", round(data.to.plot$var_accounted$Dim.1[2] * 100, 1), "% explained variation")) + 
      ylab(paste0("Axis 2: ", round(data.to.plot$var_accounted$Dim.2[2] * 100, 1), "% explained variation")) + 
      coord_cartesian(xlim = c(-max(abs(quantile(data.to.plot$ind_temp$Dim.1, c(0.001, 0.999)))),
                               max(abs(quantile(data.to.plot$ind_temp$Dim.1, c(0.001, 0.999))))),
                      ylim = c(-max(abs(quantile(data.to.plot$ind_temp$Dim.2, c(0.001, 0.999)))),
                               max(abs(quantile(data.to.plot$ind_temp$Dim.2, c(0.001, 0.999)))))) +
      theme_bw() +
      theme(panel.grid = element_blank())}
  
  if(nfactors.temp >= 4){
    plots[["pca34"]] = ggplot() +
      geom_hex(data = data.to.plot$ind_temp, aes(x = Dim.3, y = Dim.4), 
               bins = 100) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_segment(data = data.to.plot$var_temp, 
                   aes(x = 0, y = 0, xend = Dim.3 * scale.segment.dim3, yend = Dim.4 * scale.segment.dim4),
                   arrow = arrow(length = unit(0.8, "cm"))) +
      geom_label(data = data.to.plot$var_temp, aes(x = Dim.3 * scale.segment.dim3 * 1.2, y = Dim.4 * scale.segment.dim4 * 1.2, label = trait), 
                 alpha = 0.8) +
      scale_fill_gradientn(colours=c("#E8E8E8","#808080"),guide = "none", na.value = NA) +
      xlab(paste0("Axis 3: ", round(data.to.plot$var_accounted$Dim.3[2] * 100, 1), "% explained variation")) + 
      ylab(paste0("Axis 4: ", round(data.to.plot$var_accounted$Dim.4[2] * 100, 1), "% explained variation")) + 
      coord_cartesian(xlim = c(-max(abs(quantile(data.to.plot$ind_temp$Dim.3, c(0.001, 0.999)))),
                               max(abs(quantile(data.to.plot$ind_temp$Dim.3, c(0.001, 0.999))))),
                      ylim = c(-max(abs(quantile(data.to.plot$ind_temp$Dim.4, c(0.001, 0.999)))),
                               max(abs(quantile(data.to.plot$ind_temp$Dim.4, c(0.001, 0.999)))))) +
      theme_bw() +
      theme(panel.grid = element_blank())}
  
  if(nfactors.temp >= 6){
    plots[["pca56"]] = ggplot() +
      geom_hex(data = data.to.plot$ind_temp, aes(x = Dim.5, y = Dim.6), 
               bins = 100) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_segment(data = data.to.plot$var_temp, 
                   aes(x = 0, y = 0, xend = Dim.5 * scale.segment.dim5, yend = Dim.6 * scale.segment.dim6),
                   arrow = arrow(length = unit(0.8, "cm"))) +
      geom_label(data = data.to.plot$var_temp, aes(x = Dim.5 * scale.segment.dim5 * 1.2, y = Dim.6 * scale.segment.dim6 * 1.2, label = trait), 
                 alpha = 0.8) +
      scale_fill_gradientn(colours=c("#E8E8E8","#808080"),guide = "none", na.value = NA) +
      xlab(paste0("Axis 5: ", round(data.to.plot$var_accounted$Dim.5[2] * 100, 1), "% explained variation")) + 
      ylab(paste0("Axis 6: ", round(data.to.plot$var_accounted$Dim.6[2] * 100, 1), "% explained variation")) + 
      coord_cartesian(xlim = c(-max(abs(quantile(data.to.plot$ind_temp$Dim.5, c(0.001, 0.999)))),
                               max(abs(quantile(data.to.plot$ind_temp$Dim.5, c(0.001, 0.999))))),
                      ylim = c(-max(abs(quantile(data.to.plot$ind_temp$Dim.6, c(0.001, 0.999)))),
                               max(abs(quantile(data.to.plot$ind_temp$Dim.6, c(0.001, 0.999)))))) +
      theme_bw() +
      theme(panel.grid = element_blank())}
  
  if(nfactors.temp >= 8){
    plots[["pca78"]] = ggplot() +
      geom_hex(data = data.to.plot$ind_temp, aes(x = Dim.7, y = Dim.8), 
               bins = 100) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_segment(data = data.to.plot$var_temp, 
                   aes(x = 0, y = 0, xend = Dim.7 * scale.segment.dim7, yend = Dim.8 * scale.segment.dim8),
                   arrow = arrow(length = unit(0.8, "cm"))) +
      geom_label(data = data.to.plot$var_temp, aes(x = Dim.7 * scale.segment.dim7 * 1.2, y = Dim.8 * scale.segment.dim8 * 1.2, label = trait), 
                 alpha = 0.8) +
      scale_fill_gradientn(colours=c("#E8E8E8","#808080"),guide = "none", na.value = NA) +
      xlab(paste0("Axis 7: ", round(data.to.plot$var_accounted$Dim.7[2] * 100, 1), "% explained variation")) + 
      ylab(paste0("Axis 8: ", round(data.to.plot$var_accounted$Dim.8[2] * 100, 1), "% explained variation")) + 
      coord_cartesian(xlim = c(-max(abs(quantile(data.to.plot$ind_temp$Dim.7, c(0.001, 0.999)))),
                               max(abs(quantile(data.to.plot$ind_temp$Dim.7, c(0.001, 0.999))))),
                      ylim = c(-max(abs(quantile(data.to.plot$ind_temp$Dim.8, c(0.001, 0.999)))),
                               max(abs(quantile(data.to.plot$ind_temp$Dim.8, c(0.001, 0.999)))))) +
      theme_bw() +
      theme(panel.grid = element_blank())}
  
  results = list(pca_unrotated = pca.unrotated.temp,
                 pca_rotated = pca.rotated.temp,
                 data_to_plot = data.to.plot,
                 plots = plots)
  
  if(print.plot == T){
    print(results$plots$pca12)
  }
  return(results)
}

###################################
make.a.nice.pca.plot.test = function(
    data.temp, nfactors.temp = 8, rotate.temp = "varimax", print.plot = F, recode.trait.names, switch.axes.temp = NULL){
  
  pca.unrotated.temp = psych::principal(r = data.temp, nfactors = nfactors.temp, rotate = "none")
  pca.rotated.temp = psych::principal(r = data.temp, nfactors = nfactors.temp, rotate = rotate.temp)
  
  data.to.plot = list(ind_temp = as_tibble(pca.rotated.temp$scores) %>%
                        setNames(paste("Dim.", 1:nfactors.temp, sep = "")),
                      
                      var_temp = data.frame(unclass(pca.rotated.temp$loadings)) %>%
                        rownames_to_column(var = "trait") %>% 
                        as_tibble() %>%
                        setNames(c("trait", paste("Dim.", 1:nfactors.temp, sep = ""))),
                      var_accounted = data.frame(pca.rotated.temp$Vaccounted) %>%
                        setNames(paste("Dim.", 1:nfactors.temp, sep = "")))
  
  # rotate where necessary
  if(!is_empty(switch.axes.temp)){
    for(i in 1:length(switch.axes.temp)){
      data.to.plot$ind_temp[,switch.axes.temp] = - data.to.plot$ind_temp[,switch.axes.temp]
      data.to.plot$var_temp[,switch.axes.temp] = - data.to.plot$var_temp[,switch.axes.temp]}}
  
  # shorten trait names
  data.to.plot$var_temp$trait = recode.trait.names(data.to.plot$var_temp$trait)
  
  if(nfactors.temp >= 1){
    scale.segment.dim1 = 0.8 * max(abs(quantile(data.to.plot$ind_temp$Dim.1, c(0.001, 0.999)))) / 
      max(abs(data.to.plot$var_temp$Dim.1))}
  if(nfactors.temp >= 2){
    scale.segment.dim2 = 0.8 * max(abs(quantile(data.to.plot$ind_temp$Dim.2, c(0.001, 0.999)))) / 
      max(abs(data.to.plot$var_temp$Dim.2))}
  if(nfactors.temp >= 3){
    scale.segment.dim3 = 0.8 * max(abs(quantile(data.to.plot$ind_temp$Dim.3, c(0.001, 0.999)))) / 
      max(abs(data.to.plot$var_temp$Dim.3))}
  if(nfactors.temp >= 4){
    scale.segment.dim4 = 0.8 * max(abs(quantile(data.to.plot$ind_temp$Dim.4, c(0.001, 0.999)))) / 
      max(abs(data.to.plot$var_temp$Dim.4))}
  if(nfactors.temp >= 5){
    scale.segment.dim5 = 0.8 * max(abs(quantile(data.to.plot$ind_temp$Dim.5, c(0.001, 0.999)))) / 
      max(abs(data.to.plot$var_temp$Dim.5))}
  if(nfactors.temp >= 6){
    scale.segment.dim6 = 0.8 * max(abs(quantile(data.to.plot$ind_temp$Dim.6, c(0.001, 0.999)))) / 
      max(abs(data.to.plot$var_temp$Dim.6))}
  if(nfactors.temp >= 7){
    scale.segment.dim7 = 0.8 * max(abs(quantile(data.to.plot$ind_temp$Dim.7, c(0.001, 0.999)))) / 
      max(abs(data.to.plot$var_temp$Dim.7))}
  if(nfactors.temp >= 8){
    scale.segment.dim8 = 0.8 * max(abs(quantile(data.to.plot$ind_temp$Dim.8, c(0.001, 0.999)))) / 
      max(abs(data.to.plot$var_temp$Dim.8))}
  
  # remove "CWM_"
  data.to.plot$var_temp$trait = gsub("CWM_", "", data.to.plot$var_temp$trait)
  
  # remove labels with < 0.4 on two axes
  data.to.plot$var_temp_complete = data.to.plot$var_temp
  
  data.to.plot$var_temp = data.to.plot$var_temp %>% 
    mutate(Dim.1 = ifelse(abs(Dim.1) < 0.4 & abs(Dim.2) < 0.4, NA, Dim.1),
           Dim.2 = ifelse(abs(Dim.1) < 0.4 & abs(Dim.2) < 0.4, NA, Dim.2),
           Dim.3 = ifelse(abs(Dim.3) < 0.4 & abs(Dim.4) < 0.4, NA, Dim.3),
           Dim.4 = ifelse(abs(Dim.3) < 0.4 & abs(Dim.4) < 0.4, NA, Dim.4),
           Dim.5 = ifelse(abs(Dim.5) < 0.4 & abs(Dim.6) < 0.4, NA, Dim.5),
           Dim.6 = ifelse(abs(Dim.5) < 0.4 & abs(Dim.6) < 0.4, NA, Dim.6),
           Dim.7 = ifelse(abs(Dim.7) < 0.4 & abs(Dim.8) < 0.4, NA, Dim.7),
           Dim.8 = ifelse(abs(Dim.7) < 0.4 & abs(Dim.8) < 0.4, NA, Dim.8))
  
  plots = list()
  
  if(nfactors.temp >= 2){
    plots[["pca12"]] = ggplot() +
      geom_hex(data = data.to.plot$ind_temp, aes(x = Dim.1, y = Dim.2), 
               bins = 100) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_segment(data = data.to.plot$var_temp_complete, 
                   aes(x = 0, y = 0, xend = Dim.1 * scale.segment.dim1*1.1, yend = Dim.2 * scale.segment.dim2*1.1),
                   arrow = arrow(length = unit(0.8, "cm"))) +
      geom_label(data = data.to.plot$var_temp, aes(x = Dim.1 * scale.segment.dim1 * 1.2, y = Dim.2 * scale.segment.dim2 * 1.2, label = trait), 
                 alpha = 0.8) +
      scale_fill_gradientn(colours=c("#E8E8E8","#808080"),guide = "none", na.value = NA) +
      xlab(paste0("Axis 1: ", round(data.to.plot$var_accounted$Dim.1[2] * 100, 1), "% explained variation")) + 
      ylab(paste0("Axis 2: ", round(data.to.plot$var_accounted$Dim.2[2] * 100, 1), "% explained variation")) + 
      coord_cartesian(xlim = c(-max(abs(quantile(data.to.plot$ind_temp$Dim.1, c(0.001, 0.999)))),
                               max(abs(quantile(data.to.plot$ind_temp$Dim.1, c(0.001, 0.999))))),
                      ylim = c(-max(abs(quantile(data.to.plot$ind_temp$Dim.2, c(0.001, 0.999)))),
                               max(abs(quantile(data.to.plot$ind_temp$Dim.2, c(0.001, 0.999)))))) +
      theme_bw() +
      theme(panel.grid = element_blank())}
  
  if(nfactors.temp >= 4){
    plots[["pca34"]] = ggplot() +
      geom_hex(data = data.to.plot$ind_temp, aes(x = Dim.3, y = Dim.4), 
               bins = 100) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_segment(data = data.to.plot$var_temp_complete, 
                   aes(x = 0, y = 0, xend = Dim.3 * scale.segment.dim3, yend = Dim.4 * scale.segment.dim4),
                   arrow = arrow(length = unit(0.8, "cm"))) +
      geom_label(data = data.to.plot$var_temp, aes(x = Dim.3 * scale.segment.dim3 * 1.2, y = Dim.4 * scale.segment.dim4 * 1.2, label = trait), 
                 alpha = 0.8) +
      scale_fill_gradientn(colours=c("#E8E8E8","#808080"),guide = "none", na.value = NA) +
      xlab(paste0("Axis 3: ", round(data.to.plot$var_accounted$Dim.3[2] * 100, 1), "% explained variation")) + 
      ylab(paste0("Axis 4: ", round(data.to.plot$var_accounted$Dim.4[2] * 100, 1), "% explained variation")) + 
      coord_cartesian(xlim = c(-max(abs(quantile(data.to.plot$ind_temp$Dim.3, c(0.001, 0.999)))),
                               max(abs(quantile(data.to.plot$ind_temp$Dim.3, c(0.001, 0.999))))),
                      ylim = c(-max(abs(quantile(data.to.plot$ind_temp$Dim.4, c(0.001, 0.999)))),
                               max(abs(quantile(data.to.plot$ind_temp$Dim.4, c(0.001, 0.999)))))) +
      theme_bw() +
      theme(panel.grid = element_blank())}
  
  if(nfactors.temp >= 6){
    plots[["pca56"]] = ggplot() +
      geom_hex(data = data.to.plot$ind_temp, aes(x = Dim.5, y = Dim.6), 
               bins = 100) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_segment(data = data.to.plot$var_temp_complete, 
                   aes(x = 0, y = 0, xend = Dim.5 * scale.segment.dim5, yend = Dim.6 * scale.segment.dim6),
                   arrow = arrow(length = unit(0.8, "cm"))) +
      geom_label(data = data.to.plot$var_temp, aes(x = Dim.5 * scale.segment.dim5 * 1.2, y = Dim.6 * scale.segment.dim6 * 1.2, label = trait), 
                 alpha = 0.8) +
      scale_fill_gradientn(colours=c("#E8E8E8","#808080"),guide = "none", na.value = NA) +
      xlab(paste0("Axis 5: ", round(data.to.plot$var_accounted$Dim.5[2] * 100, 1), "% explained variation")) + 
      ylab(paste0("Axis 6: ", round(data.to.plot$var_accounted$Dim.6[2] * 100, 1), "% explained variation")) + 
      coord_cartesian(xlim = c(-max(abs(quantile(data.to.plot$ind_temp$Dim.5, c(0.001, 0.999)))),
                               max(abs(quantile(data.to.plot$ind_temp$Dim.5, c(0.001, 0.999))))),
                      ylim = c(-max(abs(quantile(data.to.plot$ind_temp$Dim.6, c(0.001, 0.999)))),
                               max(abs(quantile(data.to.plot$ind_temp$Dim.6, c(0.001, 0.999)))))) +
      theme_bw() +
      theme(panel.grid = element_blank())}
  
  if(nfactors.temp >= 8){
    plots[["pca78"]] = ggplot() +
      geom_hex(data = data.to.plot$ind_temp, aes(x = Dim.7, y = Dim.8), 
               bins = 100) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_segment(data = data.to.plot$var_temp_complete, 
                   aes(x = 0, y = 0, xend = Dim.7 * scale.segment.dim7, yend = Dim.8 * scale.segment.dim8),
                   arrow = arrow(length = unit(0.8, "cm"))) +
      geom_label(data = data.to.plot$var_temp, aes(x = Dim.7 * scale.segment.dim7 * 1.2, y = Dim.8 * scale.segment.dim8 * 1.2, label = trait), 
                 alpha = 0.8) +
      scale_fill_gradientn(colours=c("#E8E8E8","#808080"),guide = "none", na.value = NA) +
      xlab(paste0("Axis 7: ", round(data.to.plot$var_accounted$Dim.7[2] * 100, 1), "% explained variation")) + 
      ylab(paste0("Axis 8: ", round(data.to.plot$var_accounted$Dim.8[2] * 100, 1), "% explained variation")) + 
      coord_cartesian(xlim = c(-max(abs(quantile(data.to.plot$ind_temp$Dim.7, c(0.001, 0.999)))),
                               max(abs(quantile(data.to.plot$ind_temp$Dim.7, c(0.001, 0.999))))),
                      ylim = c(-max(abs(quantile(data.to.plot$ind_temp$Dim.8, c(0.001, 0.999)))),
                               max(abs(quantile(data.to.plot$ind_temp$Dim.8, c(0.001, 0.999)))))) +
      theme_bw() +
      theme(panel.grid = element_blank())}
  
  results = list(pca_unrotated = pca.unrotated.temp,
                 pca_rotated = pca.rotated.temp,
                 data_to_plot = data.to.plot,
                 plots = plots)
  
  if(print.plot == T){
    print(results$plots$pca12)
  }
  return(results)
}

#########################
replace.strip.colors  = function(plot.temp, strip.position, color.gradient.temp){
  
  g <- ggplot_gtable(ggplot_build(plot.temp))
  stripr <- which(grepl('strip-l', g$layout$name))
  fills <- color.gradient.temp
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  return(as_ggplot(g))
}
