# Teresita M. Porter, April 1, 2021

library(ggplot2)
library(ggmap)

# from Table S1
sites = data.frame(site=c("Clair15", "Laurel7", "Beaver18", "Clair12"),
                   lat=c(43.462901, 43.470727, 43.492001, 43.465409), 
                   lon=c(-80.58468, -80.55627, -80.60978, -80.57132), 
                   Condition=c("Good","Fair","Good","Fair"))

# Save it for later
write.csv(sites, "Sites.csv", row.names=FALSE, quote=FALSE)

base = get_map(location=c(min(sites$lon)-0.007, min(sites$lat)-0.015, 
                          max(sites$lon)+0.015, max(sites$lat)+0.015), 
               zoom=14, 
               maptype="terrain-background")

map <- ggmap(base) +
  geom_point(data=sites, aes(x=lon, y=lat, color=Condition), cex=2.5) + # plot the points
  geom_text(data=sites, aes(x=lon, y=lat, label=`site`), hjust=0.5, vjust=-0.5, size=2.5) +
  labs(x="Longitude", y="Latitude", title="Sites") + # label the axes
  theme_bw() + 
  theme(legend.position="bottom", 
        text = element_text(size=7),
        axis.text = element_text(size = rel(0.75)), 
        legend.key = element_rect(colour = "white"), 
        axis.text.x = element_text(angle=45, vjust=0.5)) # tweak the plot's appearance and legend position
map

ggsave("FigS1_Map.jpg", map, height = 4, width = 4)
# This uses stamen maps so shouldn't be any copy right issues
# Terrain Attribution to add to bottom of plot, or maybe just to legend
# Map tiles by Stamen Design, under CC BY 3.0. Data by OpenStreetMap, under ODbL.
# http://maps.stamen.com/#watercolor/12/37.7706/-122.3782
