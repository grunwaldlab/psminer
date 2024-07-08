#' Build maps
#'
#' Sample maps built using Leaflet and lat/lon coordinates in input spreadsheet
#'
#' @param sample_data
#'
#' @return Maps
#'
#' @export
library(leaflet)
make_maps <- function(sample_data){
  map_df <- data.frame(lon_coordinates = sample_data$lon_coordinates,
                       lat_coordinates = sample_data$lat_coordinates,
                       type = sample_data[,sample_data$color_by[1]])

  # Different color palette options.
  #pal <- colorFactor(palette = c('red', 'blue', 'green', 'purple', 'orange'), domain = map_df$type)
  pal <- colorFactor(palette = 'Dark2', domain = map_df$type)

  # If color_by exists, use color_by
  if(length(sample_data$color_by)>0){
    leaflet(map_df) %>% addTiles(options=providerTileOptions(noWrap = TRUE))  %>% setView(0,0, zoom = 1.4) %>%
      addCircles(lng = map_df$lon_coordinates, lat = map_df$lat_coordinates, weight=10, color = ~pal(type), label = ~as.character(map_df$type))
  } else {
    leaflet(map_df) %>% addTiles(options=providerTileOptions(noWrap = TRUE))  %>% setView(0,0, zoom = 1) %>%
      addCircles(lng = map_df$lon_coordinates, lat = map_df$lat_coordinates, weight=10)
  }
}
