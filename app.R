setwd("~/Downloads/6. Fall 2024/Skyspatial/Saudi-Arabia-Shiny 2")
library(readxl)
library(sf)
library(ggplot2)
library(rgeoboundaries)
library(rnaturalearth)
library(dplyr)
library(SpatialEpi)
library(spdep)
library(INLA)
library(sp)
library(Matrix)
library(spData)
library(akima)
library(terra)
library(raster)
library(leaflet)
library(shiny)
library(rmarkdown)
library(htmltools)



##################
# Get the map of Saudi Arabia
##################

map <- ne_states(country = "Saudi Arabia", returnclass = "sf")
map$name[5] <- "Asir" # change spellings
map <- st_set_crs(map, 4326) # certain projection
map <- map[order(map$name),] # alphabetically ordered
map <- map %>% dplyr::select(name, geometry)
maplatlong <- st_transform(map, 4326)

#################################
# Population of Saudi Arabia
#################################

# Read the Saudi Arabia population raster for the year 2017
r <- terra::rast("pop2017.tif")
# Ensure raster is in EPSG:4326 
r <- terra::project(r, "EPSG:4326")

# Aggregate raster to lower resolution
r <- terra::aggregate(r, fact=10, fun=sum)

# Sum population in regions of map
map$popraster <- terra::extract(r, map, sum, na.rm = TRUE)$pop2017

# Read the pm2.5 data as grid for 2017
data <- terra::rast("pm2.5.2017.tif")
data <- terra::aggregate(data, fact=10, fun=mean)


bb <- st_bbox(maplatlong)
data <- crop(data, bb)

# Sum pm in regions of map
map$pm2.5 <- terra::extract(data, map, mean, na.rm = TRUE)$pm2.5.2017

# Modify the data for plot
data2 <- as.data.frame(data, xy = TRUE)
names(data2) <- c("lon", "lat", "value")

dpol <- st_as_sf(data2, coords = c("lon", "lat"), crs = 4326)
dpol <- st_filter(dpol, map)
dpol$pm252012g <- as.numeric(dpol$value)

#
#


# (1) plot the population raster
r2 <- as.data.frame(r, xy = TRUE)
names(r2) <- c("lon", "lat", "value")
r2 <- st_as_sf(r2, coords = c("lon", "lat"), crs = 4326)

theme_set(theme_minimal())
boundaryregion.c <- st_union(map)
r2 <- st_filter(r2, map)
r2$populatin <- as.numeric(r2$value)
r2 <- r2[r2$populatin >= 1, ]
r2$pop <- log(r2$populatin)

# retrieve coordinates in matrix form, from an sf object
coord.r2 <- st_coordinates(r2)

# get a data frame 
DF.r2 <- data.frame(lon = coord.r2[,1], lat = coord.r2[,2], value = r2$pop)

# Interpolate to a regular grid


interp.r2 <- with(DF.r2, interp(x = lon, y = lat, z = value, 
                                xo = seq(min(lon), max(lon), length = 100), 
                                yo = seq(min(lat), max(lat), length = 100)))



# Create a raster from the interpolated result

raster.r2 <- raster(interp.r2)

# Load the Saudi Arabia boundary
boundaryregion <- geoboundaries("Saudi Arabia")
boundaryregion <- st_set_crs(boundaryregion, 4326)
boundaryregion <- st_transform(boundaryregion, crs(proj4string(raster.r2)))

# Convert boundaryregion to SpatialPolygons
boundaryregion_sp <- as(boundaryregion, "Spatial")
# Clip the raster to the Saudi Arabia boundary
clipped.raster.r2 <- mask(raster.r2, boundaryregion_sp)


range_values <- range(values(clipped.raster.r2), na.rm = TRUE)
pal <- colorNumeric("viridis", range_values, na.color = "transparent")
leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(clipped.raster.r2, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright", pal = pal, values = values(clipped.raster.r2), title = htmltools::HTML("Log(P)")) %>%
  addScaleBar(position = c("bottomleft"))





# (2) plot the PM2.5

# retrieve coordinates in matrix form, from an sf object
coord.dpol <- st_coordinates(dpol)

# get a data frame 
DF.dpol <- data.frame(lon = coord.dpol[,1], lat = coord.dpol[,2], value = dpol$pm252012g)

# Interpolate to a regular grid

interp.dpol <- with(DF.dpol, interp(x = lon, y = lat, z = value, 
                                    xo = seq(min(lon), max(lon), length = 100), 
                                    yo = seq(min(lat), max(lat), length = 100)))

# Create a raster from the interpolated result
raster.dpol <- raster(interp.dpol)

# Load the Saudi Arabia boundary
boundaryregion <- geoboundaries("Saudi Arabia")
boundaryregion <- st_set_crs(boundaryregion, 4326)
boundaryregion <- st_transform(boundaryregion, crs(proj4string(raster.dpol)))

# Convert boundaryregion to SpatialPolygons
boundaryregion_sp <- as(boundaryregion, "Spatial")
# Clip the raster to the Saudi Arabia boundary
clipped.raster.dpol <- mask(raster.dpol, boundaryregion_sp)


range_values <- range(values(clipped.raster.dpol), na.rm = TRUE)
pal <- colorNumeric("viridis", range_values, na.color = "transparent")
leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(clipped.raster.dpol, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright", bins = c(0,20,40,60,80), pal = pal, values = values(clipped.raster.dpol), title = htmltools::HTML("PM<sub>2.5</sub>")) %>%
  addScaleBar(position = c("bottomleft"))



##################
# weighted PM25 
##################


# PM25 weighted by population pm25 * pop / sum(pop)
# (PM25 at the coordinates PM25 raster) * (Population at coordinates PM25 raster)

dpol$pm25Xpop <- as.vector(dpol$pm252012g * terra::extract(r, st_coordinates(dpol)))[[1]]
# (Population at coordinates PM25 raster)

dpol$pop <- as.vector(terra::extract(r, st_coordinates(dpol)))[[1]]
# Mean air pollution within counties weighted by population
inter <- st_intersects(map, dpol)
map$pm25xpop <- sapply(inter, FUN = function(x){sum(dpol[x, ]$pm25Xpop, na.rm = TRUE)})
map$sumpop <- sapply(inter, FUN = function(x){sum(dpol[x, ]$pop, na.rm = TRUE)})
map$pm25weightedpop <- map$pm25xpop/map$sumpop

# (3) plot weighted PM25

# boundaryregion.c <- st_union(map)
# map1 <- st_as_sf(map)
#  ggplot(data = boundaryregion.c) +
#   geom_sf() +
#   geom_sf(data = map1, aes(fill = pm25weightedpop)) +
#   scale_fill_viridis_c(option = "viridis") +
#   labs(fill = expression(WPM[2.5]))

library(leaflet)
map$pm25weightedpop <- map$pm2.5 
l <- leaflet(map) %>% addTiles()
pal <- colorNumeric(palette = "viridis", domain = map$pm25weightedpop)
l %>% addPolygons(color = "grey", weight = 1, fillColor = ~pal(pm25weightedpop), fillOpacity = 0.5) %>%
  addLegend(pal = pal, values = ~pm25weightedpop, opacity = 0.5, title = htmltools::HTML("WPM<sub>2.5</sub>"), position = "bottomright")


##################
# Get the data
##################

data <- read.csv("fulldata.csv", header = T)

# Load necessary libraries
library(tidyr)
library(dplyr)

# Sample data
data <- data.frame(
  name = c("Al Bahah", "Al Hudud ash Shamaliyah", "Al Jawf", "Al Madinah", "Al Quassim", "Ar Riyad", "Ash Sharqiyah", "Asir", "Ha'il", "Jizan", "Makkah", "Najran", "Tabuk",
           "Al Bahah", "Al Hudud ash Shamaliyah", "Al Jawf", "Al Madinah", "Al Quassim", "Ar Riyad", "Ash Sharqiyah", "Asir", "Ha'il", "Jizan", "Makkah", "Najran", "Tabuk",
           "Al Bahah", "Al Hudud ash Shamaliyah", "Al Jawf", "Al Madinah", "Al Quassim", "Ar Riyad", "Ash Sharqiyah", "Asir", "Ha'il", "Jizan", "Makkah", "Najran", "Tabuk",
           "Al Bahah", "Al Hudud ash Shamaliyah", "Al Jawf", "Al Madinah", "Al Quassim", "Ar Riyad", "Ash Sharqiyah", "Asir", "Ha'il", "Jizan", "Makkah", "Najran", "Tabuk"),
  Cancer = c(777, 535, 564, 1055, 514, 12918, 6656, 2749, 996, 2602, 12423, 508, 535, 789, 10, 576, 11, 526, 12930, 6668, 2761, 1008, 2614, 12435, 520, 547, 787, 545, 574, 1065, 524, 1000, 6666, 2759, 1006, 2612, 12433, 518, 545, 791, 549, 578, 1069, 528, 12932, 6670, 2763, 100, 2616, 120, 522, 50),
  WPM2.5 = c(33.66812437, 35.68583657, 27.1143477, 36.94669693, 44.44070336, 49.93185277, 53.37439377, 31.73037286, 36.62079704, 37.73400696, 46.66303372, 41.1413486, 24.95984667,
             33.66812437, 35.68583657, 27.1143477, 36.94669693, 44.44070336, 49.93185277, 53.37439377, 31.73037286, 36.62079704, 37.73400696, 46.66303372, 41.1413486, 24.95984667,
             33.66812437, 35.68583657, 27.1143477, 36.94669693, 44.44070336, 49.93185277, 53.37439377, 31.73037286, 36.62079704, 37.73400696, 46.66303372, 41.1413486, 24.95984667,
             33.66812437, 35.68583657, 27.1143477, 36.94669693, 44.44070336, 49.93185277, 53.37439377, 31.73037286, 36.62079704, 37.73400696, 46.66303372, 41.1413486, 24.95984667),
  populatin = c(367804, 263357, 358978, 1578964, 1074007, 6227523, 3729044, 1623938, 523656, 1148297, 6555006, 414299, 648083,
                377804, 273357, 368978, 1588964, 1084007, 6237523, 3739044, 1633938, 533656, 1158297, 6565006, 424299, 658083,
                387804, 283357, 378978, 1598964, 1094007, 6247523, 3749044, 1643938, 543656, 1168297, 6575006, 434299, 668083,
                397804, 293357, 388978, 1608964, 1104007, 6257523, 3759044, 1653938, 553656, 1178297, 6585006, 444299, 678083),
  year = c(2001, 2001, 2001, 2001, 2001, 2001, 2001, 2001, 2001, 2001, 2001, 2001, 2001,
           2002, 2002, 2002, 2002, 2002, 2002, 2002, 2002, 2002, 2002, 2002, 2002, 2002,
           2003, 2003, 2003, 2003, 2003, 2003, 2003, 2003, 2003, 2003, 2003, 2003, 2003,
           2004, 2004, 2004, 2004, 2004, 2004, 2004, 2004, 2004, 2004, 2004, 2004, 2004)
)

# Pivot the data to wider format
data_wide <- data %>%
  pivot_wider(names_from = year, values_from = c(Cancer, WPM2.5, populatin), names_sep = "_")

# View the transformed data
data2 <- as.data.frame(data_wide)



d1 <- aggregate(
  x = data$Cancer,
  by = list(county = data$name, year = data$year),
  FUN = sum
)
names(d1) <- c("name", "year", "Y")



data <- data[order(
  data$name,
  data$year
), ]


n.strata <- 1

library(SpatialEpi)

E <- expected(
  population = data$populatin,
  cases = data$Cancer,
  n.strata = n.strata
)




nyears <- length(unique(data$year))
countiesE <- rep(unique(data$name),
                 each = nyears)

ncounties <- length(unique(data$name))
yearsE <- rep(unique(data$year),
              times = ncounties)

dE1 <- data.frame(name = countiesE, year = yearsE, E = E)


d1 <- merge(d1, dE1, by = c("name", "year"))

d1$SIR <- d1$Y / d1$E
head(d1)



dw1 <- reshape(d1,
               timevar = "year",
               idvar = "name",
               direction = "wide"
)




map <- merge(map, dw1, by.x = "name", by.y = "name")


map_sf <- st_as_sf(map)


library(tidyr)
map_sf <- gather(map_sf, year, SIR, paste0("SIR.", 2001:2004))

map_sf$year <- as.integer(substring(map_sf$year, 5, 8))

library(ggplot2)
ggplot(map_sf) + geom_sf(aes(fill = SIR)) +
  facet_wrap(~year, dir = "h", ncol = 4) +
  ggtitle("SIR") + theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_fill_gradient2(
    midpoint = 1, low = "blue", mid = "white", high = "red" 
  ) 




d1$idarea <- as.numeric(as.factor(d1$name))
d1$idarea1 <- d1$idarea
d1$idtime <- 1 + d1$year - min(d1$year)

formula <- Y ~  data$WPM2.5 +  offset(log(E))+f(idarea, model = "bym", graph = g) +
  f(idarea1, idtime, model = "iid") + idtime

library(INLA)
library(spdep)
nb <- poly2nb(map)
head(nb)

nb2INLA("map.adj", nb)
g <- inla.read.graph(filename = "map.adj")


res <- inla(formula,
            family = "poisson", data = d1, E = E,
            control.predictor = list(compute = TRUE), #compute the posteriors of the predictions
            control.compute = list(return.marginals.predictor = TRUE), verbose = T)


summary(res)


#########################
# Table of the estimate
#########################

library(pander)


pander(round(res$summary.fixed[,c(1,3,5)],2))



d1$RR <- res$summary.fitted.values[, "mean"]

d1$exc <- sapply(res$marginals.fitted.values,
                 FUN = function(marg){1 - inla.pmarginal(q = 1.1, marginal = marg)})


map_sf <- merge(
  map_sf, d1,
  by.x = c("name", "year"),
  by.y = c("name", "year")
)




newdata <- map_sf[, c("name", "year", "RR", "exc")]

library(tidyr)
library(dplyr)

# Assuming 'data_map_sf' is your sf object
# First, ensure that 'data_map_sf' only contains the relevant columns if there are extras
library(tidyr)
library(dplyr)

# Assuming 'data_map_sf' is your sf object
# Selecting necessary columns if there are extras
data_map_sf <- newdata %>%
  dplyr::select(name, year, RR, exc, geometry)

# Pivoting to wide format to separate both RR and exc
data_wide <- data_map_sf %>%
  pivot_wider(
    names_from = year,
    values_from = c(RR, exc),
    names_sep = "_",             # Separator between the name components
    names_prefix = "Year_"       # Prefix to indicate the year for clarity
  )


m <- as(data_wide, "Spatial")


# View the transformed data



l <- leaflet(m) %>% addTiles()
pal <- colorNumeric(palette = "viridis", domain = m$RR_Year_2001)
l %>% addPolygons(color = "grey", weight = 1, fillColor = ~pal(RR_Year_2001), fillOpacity = 0.5) %>%
  addLegend(pal = pal, values = ~RR_Year_2001, opacity = 0.5, title = htmltools::HTML("RR"), position = "bottomright")

l <- leaflet(m) %>% addTiles()
pal <- colorNumeric(palette = "viridis", domain = m$exc_Year_2001)
l %>% addPolygons(color = "grey", weight = 1, fillColor = ~pal(exc_Year_2001), fillOpacity = 0.5) %>%
  addLegend(pal = pal, values = ~exc_Year_2001, opacity = 0.5, title = htmltools::HTML("Cluster"), position = "bottomright")



map2 <- as(map, "Spatial")




library(shiny)
library(leaflet)
library(sp)

# Assuming data_wide is already defined in your environment
m <- as(data_wide, "Spatial")

library(shiny)
library(rmarkdown)
library(htmltools)
library(leaflet)
library(sf)
library(ggplot2)
library(dplyr)
library(webshot)
library(htmlwidgets)


library(shiny)
library(leaflet)
library(sf)
library(ggplot2)
library(dplyr)
library(webshot)
library(htmlwidgets)
library(rmarkdown)
library(htmltools)
library(terra)


ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .title-panel-custom {
        color: gray;
        font-size: 24px;
        font-weight: bold;
      }
      .sidebar-label-custom {
        margin-top: 20px;
      }
      .sidebar-content {
        display: flex;
        flex-direction: column;
        height: 100%;
      }
      .select-input-container {
        flex-grow: 1;
      }
      .result-container {
        margin-top: 20px;
        padding: 10px;
        border: 1px solid #ddd;
        background-color: #f9f9f9;
        border-radius: 5px;
      }
      .leaflet-container {
        margin: 0 auto; /* Centers the map and adds space around */
      }
      .map-container {
        margin-bottom: 15px; /* Adds space between maps */
      }
      .pm25-container {
        margin-bottom: 30px; /* Increased space between satellite plot and other content */
      }
    "))
  ),
  titlePanel(
    title = div("Assess the impact of pollution on cancer", class = "title-panel-custom")
  ),
  sidebarLayout(
    sidebarPanel(
      tags$div(class = "sidebar-content",
               tags$div(class = "select-input-container",
                       # tags$p("Please upload data to visualize its spatial distribution on the map."),
                        selectInput("selectData", "Select satellite data:", choices = c("", "NO3", "PM2.5", "CO")),
                        fileInput(inputId = "filedata", label = "Upload data. Choose csv file", accept = c(".csv")),
                        uiOutput("yearSlider"),
                        textInput("text", "Comments:", "")
               ),
               tags$div(id = "result-container",
                        uiOutput("resultsTitle"),
                        textOutput("report"),
                        downloadButton("downloadPDF", "Download PDF")
               )
      )
    ),
    mainPanel(
      uiOutput("satellitePlotOutput", class = "pm25-container"),
      uiOutput("mapOutput2"),
      uiOutput("mapOutput1"),
      # tableOutput("inlaResults"),
      uiOutput("markdown")
    )
  )
)

server <- function(input, output, session) {
  
  # Reactive expression to check for selected data
  selectedData <- reactive({
    input$selectData
  })
  
  output$satellitePlotOutput <- renderUI({
    if (selectedData() != "") {
      # Make the PM2.5 plot larger
      leafletOutput("satelliteMap", height = ifelse(selectedData() == "PM2.5", "238px", "155px"), width = ifelse(selectedData() == "PM2.5", "100%", "80%"))
    }
  })
  
  output$satelliteMap <- renderLeaflet({
    req(input$selectData)
    
    if (input$selectData == "PM2.5") {
      # Assuming clipped.raster.dpol is already created and available for PM2.5
      range_values <- range(values(clipped.raster.dpol), na.rm = TRUE)
      pal <- colorNumeric("viridis", range_values, na.color = "transparent")
      
      leaflet() %>%
        addProviderTiles(providers$CartoDB.Positron) %>%
        addRasterImage(clipped.raster.dpol, colors = pal, opacity = 0.5) %>%
        addLegend("bottomright",bins = c(0,15,30,45,60,75,90), pal = pal, values = values(clipped.raster.dpol), title = htmltools::HTML("PM<sub>2.5</sub>")) %>%
        addScaleBar(position = c("bottomleft"))
      
    } else if (input$selectData == "NO3") {
      # Assuming clipped.raster.no3 is already created and available for NO3
      range_values <- range(values(clipped.raster.no3), na.rm = TRUE)
      pal <- colorNumeric("viridis", range_values, na.color = "transparent")
      
      leaflet() %>%
        addProviderTiles(providers$CartoDB.Positron) %>%
        addRasterImage(clipped.raster.no3, colors = pal, opacity = 0.5) %>%
        addLegend("bottomright", pal = pal, values = values(clipped.raster.no3), title = htmltools::HTML("NO<sub>3</sub>")) %>%
        addScaleBar(position = c("bottomleft"))
    }
  })
  
  data <- reactive({
    req(input$filedata)  # Ensure the file input is available
    read.csv(input$filedata$datapath)
  })
  
  observeEvent(input$filedata, {
    output$yearSlider <- renderUI({
      sliderInput("year", "Select Year:", min = 2001, max = 2004, value = 2001, step = 1, sep = "")
    })
    output$resultsTitle <- renderUI({
      tags$h4("Results")
    })
    output$mapOutput1 <- renderUI({
      div(leafletOutput("mymap1", height = "238px", width = "100%"), class = "map-container")  # Adjusted size
    })
    output$mapOutput2 <- renderUI({
      div(leafletOutput("mymap2", height = "238px", width = "100%"), class = "map-container")  # Adjusted size
    })
  })
  
  output$mymap1 <- renderLeaflet({
    req(input$year)
    map2 <- m
    year <- input$year
    exc_col <- paste0("exc_Year_", year)
    pal <- colorNumeric(palette = "viridis", domain = map2[[exc_col]])
    
    leaflet(map2) %>%
      addTiles() %>%
      addPolygons(
        color = "grey",
        weight = 1,
        fillColor = ~pal(map2[[exc_col]]),
        fillOpacity = 0.5,
        highlight = highlightOptions(
          weight = 3,
          color = "#666",
          fillOpacity = 0.7,
          bringToFront = TRUE),
        label = ~paste(name, "Cluster:", map2[[exc_col]])
      ) %>%
      addLegend(
        pal = pal,
        values = ~map2[[exc_col]],
        opacity = 0.5,
        title = htmltools::HTML("Cluster"),
        position = "bottomright"
      )
  })
  
  output$mymap2 <- renderLeaflet({
    req(input$year)
    map2 <- m
    year <- input$year
    rr_col <- paste0("RR_Year_", year)
    pal <- colorNumeric(palette = "viridis", domain = map2[[rr_col]])
    
    leaflet(map2) %>%
      addTiles() %>%
      addPolygons(
        color = "grey",
        weight = 1,
        fillColor = ~pal(map2[[rr_col]]),
        fillOpacity = 0.5,
        highlight = highlightOptions(
          weight = 3,
          color = "#666",
          fillOpacity = 0.7,
          bringToFront = TRUE),
        label = ~paste(name, "RR:", map2[[rr_col]])
      ) %>%
      addLegend(
        pal = pal,
        values = ~map2[[rr_col]],
        opacity = 0.5,
        title = htmltools::HTML("RR"),
        position = "bottomright"
      )
  })
  
  output$report <- renderText({
    req(input$year)
    year <- input$year
    rr_col <- paste0("RR_Year_", year)
    highest_risk_region <- m$name[which.max(m[[rr_col]])]
    highest_risk_value <- max(m[[rr_col]], na.rm = TRUE)
    
    paste("If the pullution increase one unit, then the risk increase 64%, highlighting the significant impact of air pollution on public health.
          The highest relative risk (RR) in year", year, "is in", highest_risk_region, "with a value of", highest_risk_value)
  })
  output$inlaResults <- renderTable({
    req(input$year)
    res$summary.fixed
  })
  output$markdown <- renderUI({
    HTML(markdown::markdownToHTML(text = input$text, fragment.only = TRUE))
  })
  
  output$downloadPDF <- downloadHandler(
    filename = function() {
      paste("report", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      tempReport <- file.path(tempdir(), "report.Rmd")
      tempHtml <- file.path(tempdir(), "report.html")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      params <- list(text = input$text,
                     results = renderText({output$report}),
                     year = input$year
                     # inlaResults = res$summary.fixed
      )
      # Render HTML
      rmarkdown::render(tempReport, output_file = tempHtml,
                        params = params,
                        envir = new.env(parent = globalenv()))
      
      # Convert HTML to PDF
      pagedown::chrome_print(tempHtml, output = file)
    }
  )
}

shinyApp(ui = ui, server = server)
