library(tidyverse)
library(EchoNet2Fish)

sv <- readSVTS("./data/Sv/", oldname = NULL, newname = NULL, elimMiss = NULL,
               datevars = c("Date_S", "Date_E", "Date_M"), addyear = TRUE,
               tidy = TRUE)

sv38 <- readSVTS("./data/Sv38/", oldname = NULL, newname = NULL, elimMiss = NULL,
               datevars = c("Date_S", "Date_E", "Date_M"), addyear = TRUE,
               tidy = TRUE)%>%
  mutate(Frequency = 38)

sv38_120 <- rbind(sv, sv38)%>%
  select(Region_name, Interval, Layer, Layer_depth_min, Layer_depth_max,
         Lat_M, Lon_M, Exclude_below_line_depth_mean, Sv_mean, PRC_ABC,
        Frequency,  PRC_NASC)



TS <- readSVTS("./data/TS/", oldname = NULL, newname = NULL, elimMiss = NULL,
         datevars = c("Date_S", "Date_E", "Date_M"), addyear = TRUE,
         tidy = TRUE)%>%
  mutate(Frequency = 120)

TS38 <- readSVTS("./data/TS38/", oldname = NULL, newname = NULL, elimMiss = NULL,
               datevars = c("Date_S", "Date_E", "Date_M"), addyear = TRUE,
               tidy = TRUE)%>%
  mutate(Frequency = 38)

TS38_120 <- rbind(TS, TS38)
TS38_120$sigma <- sigmaAvg(TS38_120, TSrange=c(-60, -26))
lowerTS <- grep("X.60", colnames(TS38_120))
upperTS<- grep("X.26", colnames(TS38_120))
TS38_120$numtarg <- rowSums(TS38_120[col1:col2])
TS38_120a <- TS38_120 %>%
  select(Region_name, Interval, Layer, Frequency, sigma, numtarg)




min(TS38_120a$numtarg)
svts <- merge(sv38_120, TS38_120a, 
              by = c("Region_name", "Interval", "Layer", "Frequency"))

svts$nph <- 10000*(svts$PRC_ABC/svts$sigma)
svts$TS <- sigma2TS(svts$sigma)
svts_wide <- pivot_wider(svts, id_cols = c(Region_name, Interval, Layer,
              Layer_depth_min, Layer_depth_max), 
                         names_from = Frequency, 
                         values_from = c(nph, numtarg, PRC_ABC, PRC_NASC, sigma, TS))

svts_wide$dep <- NA
svts_wide$dep[svts_wide$Layer_depth_min<20] <- "Epilimnion"
svts_wide$dep[svts_wide$Layer_depth_min>=20 & svts_wide$Layer_depth_min<40] <- "Thermocline"
svts_wide$dep[svts_wide$Layer_depth_min>=40] <- "Hypolimnion"
svts_wide$Depth_layer <-factor(svts_wide$dep, levels=c('Epilimnion','Thermocline','Hypolimnion'))

svts_wide$nph_120[is.na(svts_wide$nph_120)]<-0
svts_wide$PRC_NASC_120[is.na(svts_wide$PRC_NASC_120)]<-0
svts_wide$numtarg_120[is.na(svts_wide$numtarg_120)]<-0

svts$Frequency.f <- factor(svts$Frequency)


densumm <- svts_wide %>% group_by(Region_name, Interval, Depth_layer)%>%
  summarise(nperha_38 = sum(nph_38), nperha_120 =sum(nph_120))
  
ggplot()+
  geom_point(data=densumm, aes(x=nperha_120, y = nperha_38))+
  geom_abline(slope=1, intercept=0)+
  facet_wrap(~Depth_layer, ncol=1, scales="free")

png("./output/target_numbers.png", units="in", width=12, height=8, res=900)

ggplot()+
  geom_histogram(data=subset(svts, numtarg<2500), aes(x = numtarg, 
  fill=Frequency.f), bins=15) + 
  theme(plot.subtitle = element_text(vjust = 1), 
    plot.caption = element_text(vjust = 1), 
    axis.line = element_line(linetype = "solid"), 
    panel.grid.major = element_line(linetype = "blank"), 
    panel.grid.minor = element_line(linetype = "blank"), 
    axis.title = element_text(size = 24), 
    axis.text = element_text(size = 24), 
    axis.text.x = element_text(size = 22), 
    axis.text.y = element_text(size = 24), 
    panel.background = element_rect(fill = NA)) +labs(x = "Number of targets", colour = "Frequency", 
    fill = "Frequency") + theme(legend.text = element_text(size = 24), 
    legend.title = element_text(size = 24))
dev.off()

png("./output/density_nph.png", units="in", width=12, height=8, res=1200)
ggplot()+
  geom_point(data=svts_wide, aes(x = nph_120, y = nph_38))+
  geom_abline(intercept = 0, slope =1)+
  facet_wrap(~Depth_layer, ncol=1, scales="free") + theme(plot.subtitle = element_text(vjust = 1), 
    plot.caption = element_text(vjust = 1), 
    axis.line = element_line(linetype = "solid"), 
    panel.grid.major = element_line(linetype = "blank"), 
    panel.grid.minor = element_line(linetype = "blank"), 
    axis.title = element_text(size = 23), 
    axis.text = element_text(size = 24), 
    axis.text.x = element_text(size = 23), 
    axis.text.y = element_text(size = 23), 
    panel.background = element_rect(fill = NA))
dev.off()

png("./output/nasc.png", units="in", width=12, height=8, res=1200)
ggplot()+
  geom_point(data=svts_wide, aes(x = PRC_NASC_120, y = PRC_NASC_38))+
  geom_abline(intercept = 0, slope =1)+
  facet_wrap(~Depth_layer, ncol=1, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
    plot.caption = element_text(vjust = 1), 
    axis.line = element_line(linetype = "solid"), 
    panel.grid.major = element_line(linetype = "blank"), 
    panel.grid.minor = element_line(linetype = "blank"), 
    axis.title = element_text(size = 23), 
    axis.text = element_text(size = 24), 
    axis.text.x = element_text(size = 24), 
    axis.text.y = element_text(size = 24), 
    panel.background = element_rect(fill = NA))
dev.off()

png("./output/TS.png", units="in", width=12, height=12, res=1200)
ggplot()+
  geom_point(data=svts_wide, aes(x = TS_120, y = TS_38))+
  geom_abline(intercept = 0, slope =1)+
  facet_wrap(~Depth_layer, ncol=1) + 
  theme(plot.subtitle = element_text(vjust = 1), 
    plot.caption = element_text(vjust = 1), 
    axis.line = element_line(linetype = "solid"), 
    panel.grid.major = element_line(linetype = "blank"), 
    panel.grid.minor = element_line(linetype = "blank"), 
    axis.title = element_text(size = 23), 
    axis.text = element_text(size = 24), 
    axis.text.x = element_text(size = 24), 
    axis.text.y = element_text(size = 24), 
    panel.background = element_rect(fill = NA))
dev.off()


