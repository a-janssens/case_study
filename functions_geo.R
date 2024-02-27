# Get postal code to nos_code mapping
get_post_to_nis_mapping <- function(){
  
  nis <- "Z:/Intego_Data/Spatial_Data/Other_Files/post_nis_2019.csv" 
  nis <- read.csv2(nis,sep=",") %>% 
    # Flanders
    filter(region_label_nl == "Vlaams G") %>% 
    select(postal_code, nis_code = municipality_nis_code) %>% 
    arrange(nis_code)
  
}

plot_smr <- function(plotvar,map,disease="LRTI"){
  library(RColorBrewer)
  library(classInt)
  nclr <- 7 # Numbers of colors to be used 
  plotclr <- rev(brewer.pal(nclr,"RdYlGn"))
  class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks = c(0,0.5,0.6,0.8,1.2, 1.5,2,15))
  colcode <- findColours(class,plotclr)
  plot(map,col=colcode,main=disease,cex.main=2)
  
  legend("topright",legend=names(attr(colcode,"table")),
         fill=attr(colcode,"palette"), cex = 1, bty="n", title="SIR")
  
}

get_spatial_polygon_flanders <- function(){
  
  library(rgdal)
  flanders.shp <- readOGR("/Users/u0121893/OneDrive - KU Leuven/Spatial_Data/Shape_Files/Flanders/Shapefile/",layer="Refgem")
  # plot(flanders.shp)
  
  print("Plot of Flanders")
  return(flanders.shp)
  
}


get_sf_map_flanders <- function(postal_codes=FALSE,include_bxl = FALSE){
  
  library(sf)
  
  if(postal_codes){
    
    map.sf <- st_read(dsn = "/Users/u0121893/OneDrive - KU Leuven/Spatial_Data/Shape_Files/Flanders/Shapefile",layer="Belgium_PC")
    
    if(include_bxl){
      
      bxl <- c(1070,1000,1020,1120,1130,1040,1050,
               1140,1083,1090,1081,1160,1030,1082,
               1060,1080,1210,1200,1150,1180,1190,
               1170)
      
      map.sf <- map.sf %>%  filter(POSTCODE %in% c(bxl,c(1500:3999), c(8000:9999)))
      
    }else{
      
      map.sf <- map.sf %>%  filter(POSTCODE %in% c(c(1500:3999), c(8000:9999)))
      
    }
    
    
  } else{
    
    map.sf <- st_read(dsn = "/Users/u0121893/OneDrive - KU Leuven/Spatial_Data/Shape_Files/Flanders/Shapefile",layer="Refgem")
    
  }
  
  return(map.sf)
  
}



smr_age <- function(df,age_breaks=seq(-1,105,15)){
  
  nis_data <- df %>% 
    mutate(age.cat = cut(age,breaks=age_breaks))
  
  ## Count per nis
  nis_data <- nis_data %>% 
    count(municipality_nis_code,age.cat,event)
  
  cases_i <- nis_data %>% 
    group_by(municipality_nis_code) %>% 
    mutate(y_i = sum(event==1)) %>% 
    distinct(municipality_nis_code,y_i)
  
  rates_g <- nis_data %>% 
    pivot_wider(names_from=event,names_prefix="event_", values_from = n,values_fill = 0) %>% 
    group_by(age.cat) %>% 
    mutate(r_i = sum(event_1)/sum(sum(event_1,event_0))) %>% 
    distinct(age.cat,r_i)
  
  expected_i <- nis_data %>% 
    group_by(municipality_nis_code,age.cat) %>% 
    mutate(n =  sum(n)) %>% 
    distinct(municipality_nis_code,age.cat,n) %>% 
    left_join(rates_g,by="age.cat") %>% 
    group_by(municipality_nis_code) %>% 
    mutate(e_i=sum(n*r_i)) %>% 
    distinct(municipality_nis_code,e_i)
  
  SMR <- cases_i %>% 
    left_join(expected_i,b="municipality_nis_code") %>% 
    mutate(smr_i = y_i/e_i) %>% 
    rename(nis=municipality_nis_code)
  
  return(SMR)
  
  
}


# Order map from west to eat by x coordinates of centeroids 
order_map_west_east <- function(sf_map){
  # Error message if sf_map 
  if(is.null(attr(sf_map,"sf_column"))) stop("sf_map is not an sf object.")
  
  # get xy coords 
  xy <- st_coordinates(st_centroid(sf_map))
  
  # add to sf_map
  sf_map$x <- xy[,"X"]
  sf_map$y <- xy[,"Y"]
  
  # arrange according to west-east
  sf_map_west_east <- sf_map[order(sf_map$x),]
  # add index to indicate west-east ordering 
  sf_map_west_east$west_east_index <- 1:length(xy[,1])
  # data frame 
  sf_map_we <- subset(as.data.frame(sf_map_west_east), select = -geometry)
  
  cat("Return the map with west_east_index column.\n")
  return(as_tibble(sf_map_we))
  
}
