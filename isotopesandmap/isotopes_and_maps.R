# Load Libraries
library(sf)
library(elevatr)
library(marmap)
library(rnaturalearth)
library(ggplot2)
library(raster)
library(gridExtra)
library(progress)
library(here)
library(dplyr)
library(ggrepel)
library(ggspatial)
library(ggpp)
library(ggpubr)

#Draw basemap
win  <- ne_countries(continent = 'asia',scale=10,returnclass='sf')
japan <- ne_states(country = "japan",returnclass='sf')
elevation <- get_elev_raster(locations = win, z = 4,clip = "locations",src='aws')
dem <- as.data.frame(elevation, xy = T) |> na.omit()
colnames(dem)[3] <- "elevation"
dem$elevation[which(dem$elevation < 0)] = 0

#Define datasets for plotting
S1 <- read.csv(here('data','S1.csv')) |> select(Name=Site_name,Lat=Latitude,Long=Longitude) |> unique()
S4 <- read.csv(here('data','S4.csv')) |> select(Name=SiteName,Lat=Latitude,Long=Longitude) |> unique()
S1$Type <- 'Bone'
S4$Type <- 'Pot'
SiteMaps <- rbind.data.frame(S4,S1)
SiteMaps$Num <- c(2,3,5,4,1,8,9,6,7,10,16,15,12,11,17,13,14,18,19)
SiteMaps <- SiteMaps[order(SiteMaps$Num),]
SiteMaps$Nudge_x <- c(-2,-2,1,1,1,0,-0.3,0,0,0,-1.5,-2,-2,-2,-2,-2,-2,-2,-2)
SiteMaps$Nudge_y <- c(-1,-1,0,0,0,1,1,0.25,1,1,-0.5,0,-2.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5)
SiteMaps$Num_loc <- paste(SiteMaps$Num, SiteMaps$Name)
Pots <- SiteMaps %>% filter (Type == "Pot")
Bones <- SiteMaps %>% filter (Type == "Bone")
pots.sf <- st_as_sf(Pots, coords=c('Long','Lat'),crs=4326)
bones.sf  <- st_as_sf(Bones, coords=c('Long','Lat'),crs=4326)


#Create map with site locations
map <- ggplot() +
	geom_tile(data=dem,aes(x=x,y=y,fill=elevation)) +
	scale_fill_etopo() +
	geom_sf(data=pots.sf, fill="black", shape=21, colour = "black", size=1) +
	geom_sf(data=bones.sf, fill="black", shape=21, colour = "black", size=1)+
	geom_text_repel(data = SiteMaps, aes(x = Long, y = Lat, label = Num), min.segment.length = 0, nudge_x =  SiteMaps$Nudge_x, nudge_y = SiteMaps$Nudge_y, size=2, check_overlap = T) +
	xlim(124,135) +
	ylim(30,40) +
	xlab('Longitude') +
	ylab('Latitude') +
	annotation_scale()+
	theme(axis.ticks = element_blank())+
	theme_bw(base_size = 8)+
	theme(strip.text.x = element_blank())+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
	theme(legend.position="none")


#Plot human isotope data
Data_hum <- read.csv(here('data','S4.csv')) 

points1 <- c("Bronze Age Korea" = 21, "Late Jomon" = 23, "Yayoi" = 21)
colours <- c("Bronze Age Korea" = "#E69F00", "Late Jomon" = "#0072B2", "Yayoi" = "#56B4E9")
sites <-c("Ohtomo, Yayoi" = "#0072B2", "Ohtomo, Jomon" = "#56B4E9", "Yoshinogari, Yayoi" = "#009E73","Bronze Age Korea" = "#E69F00") 


hum <-ggplot()+
	labs(y=expression(delta^{15}*N[collagen]*"(\u2030)"), x=expression(delta^{13}*C[collagen]*"(\u2030)"))+
	scale_x_continuous(limits=c(-25,-10))+
	scale_y_continuous(limits=c(2,17))+
	geom_point(data=Data_hum, aes(y=d15N_coll, x=d13C_coll, fill= Group), size=2, pch =21, colour="black")+
	scale_fill_manual(values=sites)+
	theme(axis.ticks = element_blank())+
	theme_bw(base_size = 8)+
	theme(strip.text.x = element_blank())+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
	theme(legend.position="none")+
	coord_fixed(ratio = 1)

# Combined and save figure

Fig1 <- ggarrange(map, hum,
		  labels = c("B", "C"),
		  ncol = 2, nrow = 1, widths = c(1, 1))
ggsave(here('figures','Fig1_bc.png'), width = 17.8, height = 17.8, dpi = 300, units = "cm", device='png')
ggsave(here('figures','Fig1_bc.svg'), width = 11.4, height = 11.4, dpi = 300, units = "cm", device='svg')
Fig1
