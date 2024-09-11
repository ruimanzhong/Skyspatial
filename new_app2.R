library(here)
library(sf)

rm(list = ls())
load(here("map_pop_pm2.5_2000_2019.RData"))
initial_map <- map
map <- initial_map %>% 
  filter(year == 2000) #%>% dplyr::select(geometry)

data <- initial_map %>% 
  st_drop_geometry() %>% 
  mutate(Cancer = sample(500:700, size = n(), replace = TRUE)) %>%
  rename(populatin = population, WPM2.5 = pm25weightedpop) 

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


map_sf <- gather(map_sf, year, SIR, paste0("SIR.", 2001:2004))

map_sf$year <- as.integer(substring(map_sf$year, 5, 8))


d1$idarea <- as.numeric(as.factor(d1$name))
d1$idarea1 <- d1$idarea
d1$idtime <- 1 + d1$year - min(d1$year)

formula <- Y ~  data$WPM2.5 +  offset(log(E))+f(idarea, model = "bym", graph = g) +
  f(idarea1, idtime, model = "iid") + idtime


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


map2 <- as(map, "Spatial")


m <- as(data_wide, "Spatial")




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
        addTiles() %>%
        addRasterImage(clipped.raster.dpol, colors = pal, opacity = 0.5) %>%
        addLegend("bottomright",bins = c(0,15,30,45,60,75,90), pal = pal, values = values(clipped.raster.dpol), title = htmltools::HTML("PM<sub>2.5</sub>")) %>%
        addScaleBar(position = c("bottomleft"))
      
    } else if (input$selectData == "NO3") {
      # Assuming clipped.raster.no3 is already created and available for NO3
      range_values <- range(values(clipped.raster.no3), na.rm = TRUE)
      pal <- colorNumeric("viridis", range_values, na.color = "transparent")
      
      leaflet() %>%
        addTiles() %>%
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
