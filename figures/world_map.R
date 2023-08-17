library(ggplot2)
library(sp)
library(maptools)
library(maps)
library(ggmap)
library(readxl)
library(leaflet)
library(dplyr)
library(tidyverse)
library(htmlwidgets)


#svg("Map-interactive.legend.svg")
df<-read_excel("world_map_input.xlsx")
head(df)
custom_icons <- iconList(
  "Malaspina" = makeIcon("template/triangle.png", iconWidth = 12, iconHeight = 12),
  "HOT" = makeIcon("template/diamond.png", iconWidth = 20, iconHeight = 20),
  "BATS" = makeIcon("template/icons8-circle-67.png", iconWidth = 20, iconHeight = 20)
)
#https://cdn-icons.flaticon.com/png/128/649/premium/649738.png?token=exp=1656855725~hmac=5a9bf62bb2b55a408a405d248a0ec11e
html_legend <- paste(
  "<img src='", 'https://cdn-icons-png.flaticon.com/128/481/481078.png', "'style='vertical-align:middle;width:12px;height:12px;'>", "<div style='display:inline-block;vertical-align:middle'><span style='font-style:italic'>Tara</span> Ocean</div>", "<br/>",
  "<img src='", 'template/triangle.png', "' style='vertical-align:middle;width:12px;height:12px;'>", "<div style='display:inline-block;vertical-align:middle'>Malaspina</div>", "<br/>",
  "<img src='", 'https://img.icons8.com/external-flat-icons-inmotus-design/2x/external-circle-grid-flat-icons-inmotus-design.png', "' style='vertical-align:middle;width:12px;height:12px;'>", "<div style='display:inline-block;vertical-align:middle'>BATS</div>", "<br/>",
  "<img src='", 'https://img.icons8.com/ios/2x/star.png', "' style='vertical-align:middle;width:12px;height:12px;'>", "<div style='display:inline-block;vertical-align:middle'> HOT</div>", "<br/>",
  "<img src='", 'https://cdn-icons-png.flaticon.com/128/136/136830.png', "' style='vertical-align:middle;width:12px;height:12px;'>", "<div style='display:inline-block;vertical-align:middle'>GEOTRACES</div>", "<br/>",
  "<img src='", 'https://cdn-icons-png.flaticon.com/128/399/399278.png', "' style='vertical-align:middle;width:12px;height:12px;'>", "<div style='display:inline-block;vertical-align:middle'><span style='font-style:italic'>Tara</span> Ocean Metatranscriptome</div>", "<br/>"
)

p = colorFactor(palette = c("#000000","#000080","#008B00","#FFFF00"),domain = c("bathypelagic","mesopelagic","deep chlorophyll Maximum","surface"),ordered = T)
leaflet_map =df %>% leaflet(options = leafletOptions(worldCopyJump = T, zoomControl = F, dragging = T)) %>% addProviderTiles(providers$Esri.OceanBasemap) %>%
  addMapPane("tara_mt", zIndex = 605) %>% 
  addMapPane("mala", zIndex = 610) %>% 
  addMapPane("tara", zIndex = 615) %>% 
  addMapPane("hots", zIndex = 625) %>%
  addMapPane("bats", zIndex = 630) %>%
  addMapPane("bgeo", zIndex = 635) %>%
  addCircleMarkers(data=df %>% filter(project=='Tara Ocean'), ~Longitude,~Latitude, color = ~color,radius = 4, fillOpacity = 1,stroke = FALSE,options = pathOptions(pane = "tara")) %>%
  addCircleMarkers(data=df %>% filter(project=='Tara Ocean_MT'), ~Longitude,~Latitude,radius = 6, color = ~color,fillOpacity = 1,stroke = FALSE,options = pathOptions(pane = "tara_mt")) %>%
  addRectangles(data=df %>% filter(project=='GEOTRACES'), lng1=~Longitude,lat1=~Latitude, lng2=~Longitude1,lat2=~Latitude1,color = ~color,fillOpacity = 1,stroke = FALSE,options = pathOptions(pane = "bgeo")) %>%
  #addMarkers(data=df %>% filter(project=='GEOTRACES'),~Longitude,~Latitude,icon = ~custom_icons[project], options = pathOptions(pane = "bgeo")) %>%
  addMarkers(data=df %>% filter(project=='Malaspina'), ~Longitude,~Latitude,icon = ~custom_icons[project], options = pathOptions(pane = "mala")) %>%
  addMarkers(data=df %>% filter(project=='HOT'), ~Longitude,~Latitude,icon = ~custom_icons[project], options = pathOptions(pane = "hots")) %>%
  addMarkers(data=df %>% filter(project=='BATS'), ~Longitude,~Latitude,icon = ~custom_icons[project], options = pathOptions(pane = "bats")) %>%
  #addMarkers(data=df %>% filter(project='Tara Ocean_MT'), ~Longitude,~Latitude, icon = ~custom_icons['Tara Ocean'], options = pathOptions(pane = "tara_mt")) %>%
  setView(lng = 0, lat = 20, zoom = 2) %>%
  addControl(html = html_legend, position = "bottomleft") %>%
  addLegend(position = "topright",opacity=1,pal = p, values = c("bathypelagic","mesopelagic","deep chlorophyll Maximum","surface"),title = "Depth")
  

leaflet_map
saveWidget(leaflet_map,file='Map-interactive-1.legend.html')
