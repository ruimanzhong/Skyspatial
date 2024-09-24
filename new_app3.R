library(readxl)
library(sf)
library(ggplot2)
library(rgeoboundaries)
library(rnaturalearth)
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
library(tidyr)
library(dplyr)
library(webshot)
library(htmlwidgets)
library(pander)
library(here)


# rm(list = ls())
# my_data_input <- read.csv(here("HepatitisC.csv"), header = TRUE, sep = ",")


data_analyzer_part1 <- function(my_data_input){

disease_df <- my_data_input %>% 
  mutate(year = year - 3) %>%
  rename(name = regions, disease = cases)

cov_names <- names(disease_df)[!names(disease_df) %in% c("name", "year", "disease")]

years_to_do_analysis <- disease_df$year %>% unique() %>% sort()
min_year <- min(years_to_do_analysis)
max_year <- max(years_to_do_analysis)


load(here("map_pop_pm2.5_2000_2019.RData"))
# These are the objects that are loaded:
clipped.raster.dpol.list <- clipped.raster.dpol.list
map <- map

return(list(min_year = min_year,
            max_year = max_year,
            clipped.raster.dpol.list = clipped.raster.dpol.list,
            map = map,
            disease_df = disease_df,
            cov_names = cov_names,
            years_to_do_analysis = years_to_do_analysis))
}



data_analyzer_part2 <- function(input_from_data_analyzer_part1){
for (name in names(input_from_data_analyzer_part1)) {assign(name, input_from_data_analyzer_part1[[name]])}
initial_map <- map
map <- initial_map %>% 
  filter(year == 2000) %>% dplyr::select(name, geometry) # 2000 is not special, any in 2000-2019 would work

data <- initial_map %>% 
  st_drop_geometry() %>% 
  filter(year %in% years_to_do_analysis) %>%
  merge(., disease_df,by = c("name", "year"))

d <- data


data <- data[order(
  data$name,
  data$year
), ]

d1 <- data.frame(name = rep(unique(data$name), each = length(unique(data$year))), 
                 year = rep(unique(data$year), times = length(unique(data$name))), 
                 E = expected(population = data$population, cases = data$disease, n.strata = 1)) %>%
  merge(d, ., by = c("name", "year")) %>%
  mutate(SIR = disease / E)


map <- reshape(d1, timevar = "year", idvar = "name", direction = "wide") %>%
  merge(map, ., by.x = "name", by.y = "name")

map_sf <- map %>% gather(year, SIR, paste0("SIR.", years_to_do_analysis)) %>%
  mutate(year = as.integer(substring(year, 5, 8)))


nb2INLA("map.adj", poly2nb(map))
g <- inla.read.graph(filename = "map.adj")

d1 <- d1 %>%
  mutate(
    idarea = as.numeric(as.factor(name)),
    idarea1 = idarea,
    idtime = 1 + year - min(year)
  )


if (length(cov_names) == 0) {
  formula <- disease ~ WPM2.5 + offset(log(E)) + f(idarea, model = "bym", graph = g) +
    f(idarea1, idtime, model = "iid") + idtime
} else {
  formula <- as.formula(paste("disease ~ WPM2.5 + ", paste(cov_names, collapse = " + "),
                              " + offset(log(E)) + f(idarea, model = \"bym\", graph = g) +
                             f(idarea1, idtime, model = \"iid\") + idtime"))
}


res <- inla(formula,
            family = "poisson", data = d1, E = E,
            control.predictor = list(compute = TRUE), #compute the posteriors of the predictions
            control.compute = list(return.marginals.predictor = TRUE), verbose = FALSE)


d1 <- d1 %>%
  mutate(
    RR = res$summary.fitted.values[, "mean"],
    exc = sapply(res$marginals.fitted.values, 
                 FUN = function(marg) { 1 - inla.pmarginal(q = 1.1, marginal = marg) })
  )

m <- merge(
  map_sf, d1,
  by.x = c("name", "year"),
  by.y = c("name", "year")
) %>%
  dplyr::select(name, year, RR, exc, geometry) %>%
  pivot_wider(
    names_from = year,
    values_from = c(RR, exc),
    names_sep = "_",             # Separator between the name components
    names_prefix = "Year_"       # Prefix to indicate the year for clarity
  )
return(list(res = res, m = m))
}


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
    title = div("Assess the impact of pollution on disease", class = "title-panel-custom")
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
  
  MY_INPUT_DATA <- reactive({
    req(input$filedata)  # Ensure the file input is available
    inFile <- input$filedata
    
    # Read the CSV file
    my_input_data <- read.csv(inFile$datapath, header = TRUE, sep = ",")
    return(my_input_data)
  })
  
  analysisResultPart1 <- reactive({
    data <- MY_INPUT_DATA()
    if (is.null(data)) {
      return(NULL)
    }
    return(data_analyzer_part1(data))
  })
  
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
  
  observeEvent(input$filedata, {
    output$yearSlider <- renderUI({
      result1 <- analysisResultPart1()
      min_year <- result1$min_year
      max_year <- result1$max_year
      
      sliderInput("year", "Select Year:", min = min_year, max = max_year, value = min_year, step = 1, sep = "")
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
  
  
  output$satelliteMap <- renderLeaflet({
    req(input$selectData)
    req(input$year)
    year <- input$year
    result1 <- analysisResultPart1()
    clipped.raster.dpol.list <- result1$clipped.raster.dpol.list
    if (input$selectData == "PM2.5") {
      range_values <- range(values(clipped.raster.dpol.list[[as.character(year)]]), na.rm = TRUE)
      pal <- colorNumeric("viridis", range_values, na.color = "transparent")
      
      leaflet() %>%
        addTiles() %>%
        addRasterImage(clipped.raster.dpol.list[[as.character(year)]], colors = pal, opacity = 0.5) %>%
        addLegend("bottomright",bins = c(0,15,30,45,60,75,90), pal = pal, values = values(clipped.raster.dpol.list[[as.character(year)]]), title = htmltools::HTML("PM<sub>2.5</sub>")) %>%
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
  
  analysisResultPart2 <- reactive({
    result1 <- analysisResultPart1()
    if (is.null(result1)) {
      return(NULL)
    }
    return(data_analyzer_part2(result1))
  })
  
  output$mymap1 <- renderLeaflet({
    req(input$year)
    year <- input$year
    result2 <- analysisResultPart2()
    m <- result2$m
    
    exc_col <- paste0("exc_Year_", year)
    pal <- colorNumeric(palette = "viridis", domain = m[[exc_col]])
    
    leaflet(m) %>%
      addTiles() %>%
      addPolygons(
        color = "grey",
        weight = 1,
        fillColor = ~pal(m[[exc_col]]),
        fillOpacity = 0.5,
        highlight = highlightOptions(
          weight = 3,
          color = "#666",
          fillOpacity = 0.7,
          bringToFront = TRUE),
        label = ~paste(name, "Cluster:", m[[exc_col]])
      ) %>%
      addLegend(
        pal = pal,
        values = ~m[[exc_col]],
        opacity = 0.5,
        title = htmltools::HTML("Cluster"),
        position = "bottomright"
      )
  })
  
  output$mymap2 <- renderLeaflet({
    req(input$year)
    year <- input$year
    result2 <- analysisResultPart2()
    m <- result2$m
    
    rr_col <- paste0("RR_Year_", year)
    pal <- colorNumeric(palette = "viridis", domain = m[[rr_col]])
    
    leaflet(m) %>%
      addTiles() %>%
      addPolygons(
        color = "grey",
        weight = 1,
        fillColor = ~pal(m[[rr_col]]),
        fillOpacity = 0.5,
        highlight = highlightOptions(
          weight = 3,
          color = "#666",
          fillOpacity = 0.7,
          bringToFront = TRUE),
        label = ~paste(name, "RR:", m[[rr_col]])
      ) %>%
      addLegend(
        pal = pal,
        values = ~m[[rr_col]],
        opacity = 0.5,
        title = htmltools::HTML("RR"),
        position = "bottomright"
      )
  })
  
  output$report <- renderText({
    req(input$year)
    year <- input$year
    result2 <- analysisResultPart2()
    res <- result2$res
    m <- result2$m
    
    rr_col <- paste0("RR_Year_", year)
    highest_risk_region <- m$name[which.max(m[[rr_col]])]
    highest_risk_value <- max(m[[rr_col]], na.rm = TRUE)
    percentage_increase <- paste0(round((highest_risk_value - floor(highest_risk_value))*100, 2), "%,")
    paste("If the pollution increases one unit, then the risk increases ", percentage_increase, "highlighting the significant impact of air pollution on public health.
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
