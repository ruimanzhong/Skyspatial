# Use the existing image with R, Shiny, and your packages installed
FROM skyspatial_inla_shiny

# Copy the Shiny app file (new_app2.R) to the container's Shiny server directory
COPY new_app2.R /srv/shiny-server/new_app2.R

# Expose port 3838 for the Shiny app
EXPOSE 3838

# Command to run the Shiny app when the container starts
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/new_app2.R', host = '0.0.0.0', port = 3838)"]

