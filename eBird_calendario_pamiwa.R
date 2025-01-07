# Aves del calendario ecológico Pamiwã en eBird

#### Paquetes ####
#eBird data filters - filtro a los datos de eBird
library(auk);
#seven packages in one for data management and visualization - siete paquetes en uno para manipulación y visualización de datos
library(tidyverse);
#Spatiotemporal Subsampling - submuestreo espaciotemporal
library(dggridR);
#load maps - cargar mapas
library(maps);
#composite figure - figuras compuestas
library(gridExtra);
library(cowplot);
#scales - escalas
library(scales)

### Datos de eBird ####

#Datos descargados eBird - actualizar a Diciembre 2024!!!
Vaupes_ebd <- "ebd_CO-VAU_smp_relNov-2024/ebd_CO-VAU_smp_relNov-2024.txt"
Vaupes_sed <- "ebd_CO-VAU_smp_relNov-2024/ebd_CO-VAU_smp_relNov-2024_sampling.txt"

#extension geográfica
pamiwa <- c(-71, 0.85,-69.84, 1.7) #W, S, E, N

#lista de especies
aves_pamiwa <- read_csv("aves_pamiwa.csv")

# cambiar a recientes cambios del nombre de algunas especies
aves_pamiwa$Especie[which(aves_pamiwa$Especie == "Vanellus cayanus")] = "Hoploxypterus cayanus"
aves_pamiwa$Especie[which(aves_pamiwa$Especie == "Celeus grammicus")] = "Celeus undatus"
aves_pamiwa$Especie[which(aves_pamiwa$Especie == "Turdus Lawrencii")] = "Turdus lawrencii"


aves_pamiwa <- aves_pamiwa$Especie

#Archivos temporales para filtro
f_ebd <- "temporal/ebd_Mitu.txt"
f_sed <- "temporal/sed_Mitu.txt"

#columnas a mantener
  #(reduce tamaño del archivo, dejando solo lo que nos interesa)
colsE <- c("observer_id", "sampling_event_identifier",
           "group identifier",
           "common_name", "scientific_name",
           "observation_count",
           "state_code", "locality_id", "latitude", "longitude",
           "protocol_type", "all_species_reported",
           "observation_date",
           "time_observations_started",
           "duration_minutes", "effort_distance_km",
           "number_observers")

#Correr filtro
ebd_filt <- auk_ebd(file = Vaupes_ebd,           #archivo de registros
                  file_sampling = Vaupes_sed) |> #archivo de muestreo
  auk_bbox(pamiwa) |>                            #extensión geográfica
  auk_species(aves_pamiwa) |>                    #lista de especies
  auk_protocol(c("Traveling", "Stationary")) |>  #los protocolos mas comunes
  auk_distance(distance = c(0,5)) |>             #eventos con distancias entre 0 y 5 km
  auk_duration(duration = c(0,300)) |>           #eventos con distancias entre 0 y 5 horas
    auk_complete() |>                            #solo listas completas
    auk_filter(f_ebd,
             f_sed,
             overwrite=T,
             keep = colsE)

#extraer los datos de muestreo
sed_only <- read_sampling(f_sed)

#extraer los datos con registros
ebd_only <- read_ebd(f_ebd)

#### detección diaria para cada especie ####

#aca vamos a iterar un proceso para hacer la historia de deteccion de las especies

# Primero,
  #una función que convierte el tiempo de observación a horas desde el amanecer
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}

#Lista desocupada para guardar los resultados
spp_list <- list()

for (spp in unique(ebd_only$scientific_name)){
  #filtrar por cada especie
  especie <- ebd_only |>
    filter(scientific_name == spp)

  #colapsar deteccion-no detección
  especie_zf <- auk_zerofill(especie, sed_only, collapse = T)

  #filtrar por metadatos de evento y seleccionar un solo registro por dia
  especie_eff <- especie_zf |>
    mutate(
      # convertir conteos a numeros
      observation_count = as.integer(observation_count),
      # poner la distancia en 0 para listas estacionarias
      effort_distance_km = if_else(protocol_type == "Stationary",
                                   0, effort_distance_km),
      # convertir el tiempo de iniciar la lista a valor decimal desde media noche
      time_observations_started = time_to_decimal(time_observations_started),
      # crear una columna donde se redondea la hora de registro a entero
      hour_sampling = round(time_observations_started, 0),
      # crear columnas para separar fecha en año, mes, semana y día del año
      year = year(observation_date),
      month = month(observation_date),
      week = week(observation_date),
      day_of_year = yday(observation_date)) |>
    # filtrar para mejorar la inferencia ecológica
    filter(number_observers <= 10,	   #solo listas con menos de 10 observadores
           effort_distance_km <= 4,    #asegurarse del esfuerzo por distancia es maximo el mismo que SENA
           duration_minutes <= 300,    #asegurarse del esfuerzo por duración
           hour_sampling %in% (5:11)) |>  #las mismas horas de muestreo que SENA
    # remover conteos sin números
    drop_na(observation_count) |>
    # agrupar por dia de año
    group_by(day_of_year) |>
    # seleccionar solo un evento, aleatoriamente, por dia
    sample_n(size = 1, replace = TRUE)

spp_list[[spp]] <- especie_eff #Add the species to the list
}

#convertir la lista en un data frame
aves_pamiwa_df <- bind_rows(spp_list)

#### Linea de tiempo de listas y registros ####

#línea de tiempo de listas y registros

fig6a <- aves_pamiwa_df |>
  mutate(detection = ifelse(observation_count >= 1,
                            1, observation_count)) |>
  group_by(checklist_id, year(observation_date), detection) |>
  summarise(n()) |>
  ggplot(aes(x = `year(observation_date)`,
             fill = factor(detection)))+
  geom_bar()+
  labs(x = "Año",
       y = "Numero de listas",
       tag = "(a)",
       fill = "Detección") +
  scale_x_continuous(breaks = seq(2011,2024,4))+
  scale_fill_manual(values = c("lightgray","black"))+
  theme_bw() +
  theme(legend.position = c(0.25,0.75),
        legend.background = element_rect(fill = "white"),
        legend.box.background = element_rect(colour = "black"))
fig6a

fig6b <- aves_pamiwa_df |>
  mutate(detection = ifelse(observation_count >= 1,
                            1, observation_count)) |>
  ggplot(aes(x = year(observation_date),
                    fill = factor(detection)))+
  geom_bar(position = "fill")+
  scale_x_continuous(breaks = seq(2011,2024,4))+
  scale_fill_manual(values = c("lightgray", "black"))+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),
                     labels = scales::percent(c(0,0.25,0.5,0.75,1)))+
  labs(x = "Año",
       y = "Porcentaje de registros",
       tag = "(b)") +
  theme_bw() +
  theme(legend.position = "none")
fig6b

#Cuantas especies por mes entre 2019-2024?

species_per_month <- aves_pamiwa_df |>
  filter(species_observed == TRUE,
         year %in% c(2019:2024)) |>
  group_by(year,
           month) %>%
  summarise(num_species = n_distinct(scientific_name))

aves_pamiwa_df <- aves_pamiwa_df |>
  left_join(species_per_month,
            by = c("year", "month"))

fig6c <- aves_pamiwa_df |>
  #filtrar solo las detecciones positivas
  filter(year %in% c(2019, 2022:2024)) |>
  #generar una columna nueva con la "deteccion" (1, detectado vs 0, no detectado)
  mutate(detection = ifelse(observation_count >= 1,
                            1, observation_count)) |>
  ggplot(aes(x = month,
             fill = factor(detection))) +
  geom_bar(position = "fill")+
  theme_bw()+
  theme(legend.position = "none",
        strip.background = element_blank()) +
  geom_segment(aes(x = 11.5, xend = 12.5, #Okorimi
                   y = Inf,
                   yend = Inf),
               color = "#B25628")+
  geom_segment(aes(x = 0.5, xend = 3, #Okorimi
                     y = Inf,
                     yend = Inf),
                 color = "#B25628")+
  geom_segment(aes(x = 3, xend = 5.5, #Okorimi
                   y = Inf,
                   yend = Inf),
               color = "#6093D4")+
  geom_segment(aes(x = 5.5, xend = 8, #Okorimi
                   y = Inf,
                   yend = Inf),
               color = "#335E96")+
  geom_segment(aes(x = 8, xend = 11.5, #Pamurimi
                   y = Inf,
                   yend = Inf),
               color = "#4C7C31")+
  geom_vline(xintercept = c(3,     #marzo 15
                              8,    #agosto 14
                              11.5),
              color = "darkgray") + #noviembre 30
  coord_polar(start = -0.5) +
  scale_x_continuous(breaks = seq(1,12,1),
                     labels = c("ene",
                                "feb",
                                "mar",
                                "abr",
                                "may",
                                "jun",
                                "jul",
                                "ago",
                                "sep",
                                "oct",
                                "nov",
                                "dic"))+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),
                     labels = scales::percent(c(0,0.25,0.5,0.75,1)))+
  scale_fill_manual(values = c("lightgray", "black"))+
  facet_wrap(~year, ncol = 2) +
  labs(x = "Años con mejor representación",
       y = "Porcentaje de registros por mes",
       tag = "(c)",
       fill = "Detección") +
  geom_text(aes(y = 0.75,
                label = num_species),
            color = "black")

#combinar las figuras
fig6ab <- plot_grid(fig6a, fig6b,
             ncol = 2)

Fig6 <- plot_grid(fig6ab, fig6c,
                  ncol = 1,
                  rel_heights = c(1/4, 1/2))

Fig6
#Save 5x7.5
ggsave(filename = "Calendario aves pamiwaeBird.pdf",
       plot = Fig6, units = "in", width = 5.5, height = 7.5)

### Especies ejemplo ####

ejemplo_spp <- c("Ramphastos tucanus",
                 "Brotogeris cyanoptera",
                 "Amazona kawalli",
                 "Tyrannus savana",
                 "Rhegmatorhina cristata",
                 "Rupicola rupicola")

fig7 <- aves_pamiwa_df |>
  filter(scientific_name %in% ejemplo_spp) |>
  group_by(scientific_name, month) |>
  summarise(mu_abund = max(observation_count, na.rm = T)) |> #View()
  ggplot(aes(x = month,
             y = as.integer(mu_abund),
             fill = scientific_name)) +
  geom_col()+
  facet_wrap(~factor(scientific_name,
                     levels = c("Ramphastos tucanus",
                                "Brotogeris cyanoptera",
                                "Amazona kawalli",
                                "Tyrannus savana",
                                "Rhegmatorhina cristata",
                                "Rupicola rupicola")),
             ncol = 3,
             scales = "free_y") +
  coord_polar(start = -0.5) +
  geom_segment(aes(x = 11.5, xend = 12.5, #Okorimi
                   y = Inf,
                   yend = Inf),
               color = "#B25628")+
  geom_segment(aes(x = 0.5, xend = 3, #Okorimi
                   y = Inf,
                   yend = Inf),
               color = "#B25628")+
  geom_segment(aes(x = 3, xend = 5.5, #Okorimi
                   y = Inf,
                   yend = Inf),
               color = "#6093D4")+
  geom_segment(aes(x = 5.5, xend = 8, #Okorimi
                   y = Inf,
                   yend = Inf),
               color = "#335E96")+
  geom_segment(aes(x = 8, xend = 11.5, #Pamurimi
                   y = Inf,
                   yend = Inf),
               color = "#4C7C31")+
  geom_vline(xintercept = c(3,     #marzo 15
                            8,    #agosto 14
                            11.5),
             color = "darkgray") + #noviembre 30
  scale_x_continuous(breaks = seq(1,12,1),
                     labels = c("ene",
                                "feb",
                                "mar",
                                "abr",
                                "may",
                                "jun",
                                "jul",
                                "ago",
                                "sep",
                                "oct",
                                "nov",
                                "dic"))+
  scale_fill_manual(values = c("#8b4aa6",
                               "#6fc14d",
                               "black",
                               "#824b35",
                               "#ff4d01",
                               "#afb1ae"))+
  scale_y_continuous(label = scales::label_number(accuracy = 1))+
  labs(x = "Seis especies (ejemplo)",
       y = "Máximo número de individuos por mes") +
  theme_bw()+
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "italic"))

fig7
ggsave("Anual_abundancia_6aves_pamiwa.png",
       fig7, width = 6.6, height = 4.5, units = "in")
