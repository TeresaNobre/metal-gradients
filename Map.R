#######################################################
#### Project: Heavy metal contamination gradients  ####
####          in a tropical seagrass ecosystem     ####
#### Script purpose: Visualisation of the site map ####
#######################################################

#### 1.  Spatial data ####
#### 1.1 Sites ####
sites <- read.csv("~/Documents/MScProject/Data/gradient.csv")
# load island data
site <- read.csv("~/Documents/MScProject/Data/comparison.csv") # load mainland data
sites <- sites[c(1,4,19,20,34,49,64,65,74,83,86,89,104,119,
                 131,134,149,161,163,165,167,168,169),2:5] # extract useful data
site <- site[1,2:4] # extract useful data

#### 1.2 Base maps ####
install.packages("rworldmap")
install.packages("rworldxtra")
require(rworldmap) # this package contains polygons for the world's landmasses
overview <- getMap(resolution = "high") # extract a high resolution world map


require(sf) # this package makes shapefile manipulation possible
# note that despite my effort to reduce their size through cropping, the following shapefiles
# take very long to load; allow for between five and ten minutes
Sper.crop <- c(xmin = 119.2, xmax = 119.5, ymin = -5.3, ymax = -4.8) # limits of the Spermonde Archipelago map
land <- read_sf("~/Documents/MScProject/Data/OSM/land_polygons.shp") %>% 
  st_crop(Sper.crop) # crop Indonesian landmasses shapefile down to the Spermonde Archipelago
roads <- read_sf("~/Documents/MScProject/Data/OSMextra/gis_osm_roads_free_1.shp") %>% 
  st_crop(Sper.crop) # crop Indonesian roads shapefile down to the Spermonde Archipelago
Lomp.crop <- c(xmin = 119.32, xmax = 119.34, ymin = -5.06, ymax = -5.04) # limits of Pulau Bara Lompo
build1 <- read_sf("~/Documents/MScProject/Data/OSMextra/gis_osm_buildings_a_free_1.shp") %>% 
  st_crop(Lomp.crop) # crop Indonesian buildings shapefile down to Pulau Bara Lompo
build2 <- read_sf("~/Documents/MScProject/Data/OSMextra/gis_osm_buildings_a_free_3.shp") %>% 
  st_crop(Lomp.crop) # crop Indonesian buildings shapefile down to Pulau Bara Lompo

#### 1.3 Bathymetry ####
# require(marmap) # this package contains global bathymetry data
# bathy <- getNOAA.bathy(lon1 = 119, lon2 = 120, lat1 = -6, lat2 = -4,
#                        resolution = 1) # extract high 30 arcsecond bathymetry data
# bathy <- fortify.bathy(bathy) # convert bathymetry data into dataframe
reef <- read.csv("~/Documents/MScProject/Data/reef.edge.csv") # load reef edge data

#### 2.  Visualisation ####
#### 2.1 Set theme ####
require(ggplot2) # this package makes elegant graphics manipulation possible
mytheme <- theme(plot.background = element_blank(), # remove plot background
                 panel.background = element_blank(), # remove panel background
                 panel.grid.major = element_blank(), # remove large panel grids
                 panel.grid.minor = element_blank(), # remove small panel grids
                 panel.border = element_rect(fill = NA, size = 1), # set panel border
                 axis.line = element_blank(), # remove axis lines
                 axis.title = element_blank(), # remove axis titls
                 axis.text = element_text(size = 12, colour = "black"), # set axis text size
                 axis.ticks.length = unit(.25, "cm"), # set axis tick length
                 axis.ticks = element_line(colour = "black"), # set axis tick colour
                 legend.key = element_blank(),  # remove legend key
                 legend.text = element_text(size = 12), # set legend text size
                 legend.text.align = 0, # align legend text
                 legend.title = element_text(size = 12, face = "bold"), # set legend title size
                 legend.background = element_blank(), # remove legend background
                 text = element_text(family = "Helvetica")) # set font

#### 2.2 Plot map of Sulawesi ####
Sulawesi <- ggplot() + # draw plotting environment
                geom_polygon(overview, mapping = aes(long, lat, group = group), 
                             colour = "#898b8e", fill = "#b5b8ba", size = 0.2) + # draw landmasses as polygons
                geom_rect(aes(xmin = 119.259, xmax = 119.441, ymin = -5.187, ymax = -4.913),
                          fill = NA, colour = "#000000", size = 0.5, linejoin = "mitre") + # draw rectangle
                coord_sf(xlim = c(116.35, 123.65), ylim = c(-5.64, 1.64)) + # set coordinate system and limits
                scale_x_continuous(breaks = seq(116, 124, by = 2),
                                   labels = c("116°E","118°E","120°E","122°E","124°E")) + # customise x-axis breaks and labels
                scale_y_continuous(breaks = seq(-6, 2, by = 2),
                                   labels = c("6°S","4°S","2°S","0°","2°N")) + # customise y-axis breaks and labels
                mytheme # add custom theme (defined above)
Sulawesi # plot

#### 2.3 Plot map of the Spermonde Archipelago and Kota Makassar ####
Spermonde <- ggplot() +
              geom_sf(land, mapping = aes(), colour = "#898b8e", fill = "#b5b8ba", size = 0.2) + # draw landmasses as simple features object
              geom_sf(roads, mapping = aes(), colour = "#000000", size = 0.1) + # draw roads as simple features object
              geom_point(site, mapping = aes(Longitude, Latitude), 
                         size = 2, colour = "#ef8407") + # draw site point
              geom_rect(aes(xmin = 119.3207, xmax = 119.3353, ymin = -5.0553, ymax = -5.0407),
                        fill = NA, colour = "#000000", size = 0.5, linejoin = "mitre") + # draw rectangle
              coord_sf(xlim = c(119.259, 119.441), ylim = c(-5.187, -4.913)) + # set coordinate system and limits
              scale_x_continuous(breaks = seq(119.25, 119.45, by = 0.1)) + # customise x-axis breaks
              scale_y_continuous(breaks = seq(-5.2, -4.9, by = 0.1)) + # customise y-axis breaks
              mytheme # add custom theme (defined above)
Spermonde # plot

#### 2.4 Plot map of Pulau Bara Lompo and the sampling locations ####
Lompo <- ggplot() +
            geom_polygon(reef, mapping = aes(Longitude, Latitude, group = Group),
                         colour = "#000000", fill = NA, size = 0.2) + # draw reef edges as polygons
            geom_sf(land, mapping = aes(), colour = NA, fill = "#b5b8ba") + # draw island as simple features object
            geom_sf(roads, mapping = aes(), colour = "#000000", size = 0.2) + # draw roads as simple features object
            geom_sf(build1, mapping = aes(), colour = NA, fill = "#898b8e") + # draw first set of buildings as simple features object
            geom_sf(build2, mapping = aes(), colour = NA, fill = "#898b8e") + # draw second set of buildings as simple features object
            geom_point(sites, mapping = aes(Longitude, Latitude),
                       size = 2, shape = 4, colour = "#ef8407") + # draw site points
            scale_y_continuous(breaks = seq(-5.056, 5.04, by = 0.004)) + # customise y-axis breaks
            scale_x_continuous(breaks = seq(119.32, 119.336, by = 0.004)) + # customise x-axis breaks
            coord_sf(xlim = c(119.3207, 119.3353), ylim = c(-5.0553, -5.0407)) + # set coordinate system and limits
            mytheme # add custom theme (defined above)
Lompo # plot

#### 2.5 Plot combined map ####
require(cowplot) # this package makes combined plotting of ggplots possible
combined <- plot_grid(Sulawesi, Spermonde, Lompo, labels = c("a", "b", "c"),
                      label_size = 12, ncol = 3, rel_widths = c(1, 0.777777777777778, 1)) # combine maps
combined # plot (optimal dimensions: 4 x 13 in)

#### 3.  Clean up ####
detach(package:rworldmap) # detach rworldmap package
detach(package:sf) # detach sf package
detach(package:ggplot2) # detach ggplot2 package
detach(package:cowplot) # detach cowplot package
rm(list = ls()) # clear environment
graphics.off() # clear plots
cat("\014") # clear console
