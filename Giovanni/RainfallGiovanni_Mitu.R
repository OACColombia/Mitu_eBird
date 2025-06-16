#Hobos and lists
rm(list=ls())

library(ggplot2) #pretty figures
library(dplyr) #%>% and filter()
library(cowplot) #ggdraw


library(tidyr) #gather et al
library(MASS) #glm.nb
library(broom) # for tidy model output
library(vegan) #NMDS

#################
# Time Series, Area-Averaged of Daily accumulated precipitation (combined microwave-IR) estimate
#################
#Final Run daily 0.1 deg. [GPM GPM_3IMERGDF v06] mm over 2000-06-01 - 2021-09-30,
  #Region 70.25W, 1N, 69.75W, 1.5N

rainfall<-read.csv("Giovanni/RainfallGiovanni.csv")
str(rainfall)

hist(rainfall$DailyPrecipitation, breaks = 50)

#time is as character... we have to change to Date.
#in addition, extract the mean precipitation value per mean per month

rainfall2 <- rainfall %>%
  mutate(Date = as.Date(time,format =  "%d/%m/%Y")) %>%
  mutate(Month = format(Date, "%m"),
         Week = week(Date))%>%
  mutate(Day = as.numeric(format(Date, "%d"))) %>%
  mutate(year = format(Date, "%Y")) %>%
  mutate(Year = as.numeric(year)+2000)%>%
  group_by(Year, Week)%>%
  summarise_at(vars(DailyPrecipitation),list(Precipitation = mean))

PrecMax = max(rainfall2$Precipitation)
PrecMean = mean(rainfall2$Precipitation)
PrecSD = sd(rainfall2$Precipitation)

rainfall2 <- rainfall2 |>
  mutate(Season = ifelse(Week %in% c(33:48),"Pamurimi",
                         ifelse(Week %in% c(49:53,1:10),"Ihirimi",
                                ifelse(Week %in% c(11:32),"Okorimi",
                                       "other"))))
#Temperature
tempera<-read.csv("Giovanni/TemperatureGiovanni.csv")
str(tempera)

hist(tempera$mean_Temperature, breaks = 50)

#time is as character... we have to change to Date.
#in addition, extract the mean value of Temperature per year per month

tempera2 <- tempera %>%
  mutate(Date = as.Date(time,format =  "%d/%m/%y")) %>%
  mutate(Month = format(Date, "%m"),
         Week = week(Date))%>%
  mutate(Day = as.numeric(format(Date, "%d"))) %>%
  mutate(Year = format(Date, "%Y")) %>%
  group_by(Year, Week)%>%
  summarise_at(vars(mean_Temperature),list(Temperature = mean))

TempMax = max(tempera2$Temperature)
TempMean = mean(tempera2$Temperature)
TempSD = sd(tempera2$Temperature)

tempera2 <- tempera2 |>
  mutate(Season = ifelse(Week %in% c(33:48),"Pamurimi",
                         ifelse(Week %in% c(49:53,1:10),"Ihirimi",
                                ifelse(Week %in% c(11:32),"Okorimi",
                                       "other"))))


#Climate change is occurring!
tempera2$year = as.numeric(tempera2$Year)

rainfall2$year = as.numeric(rainfall2$Year)

ggplot(rainfall2, aes(x=year,
                      y = Precipitation))+
  geom_point(aes(fill=factor(Season,
                             levels = c("Pamurimi",
                                        "Ihirimi",
                                        "Okorimi"))),
             shape=21,
             alpha = 0.6)+
  scale_fill_manual(values = c("#3B7D23","#C04F15","#4E95D9"),
                     labels = c("Pamurimi" = paste0("Pamur",
                                                    "\u0336","i",
                                                    "m",
                                                    "\u0336","i",
                                                    "\n (transition)"),
                                "Ihirimi" = paste0("\u0336","I",
                                                   "h",
                                                   "\u0336","i",
                                                   "r",
                                                   "\u0336","i",
                                                   "m",
                                                   "\u0336","i",
                                                   "\n (dry)"),
                                "Okorimi" = paste0("Okor",
                                                   "\u0336","i",
                                                   "m",
                                                   "\u0336","i",
                                                   "\n (rain)")))+
  geom_smooth(method = "gam", aes(color=factor(Season,
                                               levels = c("Pamurimi",
                                                          "Ihirimi",
                                                          "Okorimi"))))+
  scale_fill_manual(values = c("#3B7D23","#C04F15","#4E95D9"),
                    labels = c("Pamurimi" = paste0("Pamur",
                                                   "\u0336","i",
                                                   "m",
                                                   "\u0336","i",
                                                   "\n (transition)"),
                               "Ihirimi" = paste0("\u0336","I",
                                                  "h",
                                                  "\u0336","i",
                                                  "r",
                                                  "\u0336","i",
                                                  "m",
                                                  "\u0336","i",
                                                  "\n (dry)"),
                               "Okorimi" = paste0("Okor",
                                                  "\u0336","i",
                                                  "m",
                                                  "\u0336","i",
                                                  "\n (rain)")))+
  labs(x="",
       y = "Precipitation (mm)",
       fill = "Season",
       color = "Season")+
  theme_classic()+
  theme(legend.position = "none")


ggplot(tempera2, aes(x=year,
                      y = Temperature))+
  geom_point(aes(fill=factor(Season,
                             levels = c("Pamurimi",
                                        "Ihirimi",
                                        "Okorimi"))),
             shape=21,
             alpha = 0.6)+
  scale_fill_manual(values = c("#3B7D23","#C04F15","#4E95D9"),
                    labels = c("Pamurimi" = paste0("Pamur",
                                                   "\u0336","i",
                                                   "m",
                                                   "\u0336","i",
                                                   "\n (transition)"),
                               "Ihirimi" = paste0("\u0336","I",
                                                  "h",
                                                  "\u0336","i",
                                                  "r",
                                                  "\u0336","i",
                                                  "m",
                                                  "\u0336","i",
                                                  "\n (dry)"),
                               "Okorimi" = paste0("Okor",
                                                  "\u0336","i",
                                                  "m",
                                                  "\u0336","i",
                                                  "\n (rain)")))+
  geom_smooth(method = "gam", aes(color=factor(Season,
                                               levels = c("Pamurimi",
                                                          "Ihirimi",
                                                          "Okorimi"))))+
  scale_fill_manual(values = c("#3B7D23","#C04F15","#4E95D9"),
                    labels = c("Pamurimi" = paste0("Pamur",
                                                   "\u0336","i",
                                                   "m",
                                                   "\u0336","i",
                                                   "\n (transition)"),
                               "Ihirimi" = paste0("\u0336","I",
                                                  "h",
                                                  "\u0336","i",
                                                  "r",
                                                  "\u0336","i",
                                                  "m",
                                                  "\u0336","i",
                                                  "\n (dry)"),
                               "Okorimi" = paste0("Okor",
                                                  "\u0336","i",
                                                  "m",
                                                  "\u0336","i",
                                                  "\n (rain)")))+
  labs(x="",
       y = "Temperature (ºC)",
       fill = "Season",
       color = "Season")+
  theme_classic()+
  theme(legend.position = "none")

