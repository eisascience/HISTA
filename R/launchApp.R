#' Launch HISTA
#' @title Launch HISTA
#' @description Launch HISTA
#' @keywords shiny SDA Testis Single-Cell RNASeq
#' @export
#' @return Shiny application.
#' @import shiny
#' @import shinydashboard
#' @import shinyWidgets
#' @import shinyFiles
#' @import data.table
#' @rawNamespace import(dplyr, except = c(last, first, select, between, combine)) 
#' @import ggplot2
#' @import ggforce
#' @import ggthemes
#' @import ggpubr
#' @import ggnewscale
#' @import ggrepel
#' @import rclipboard
#' @import grid
#' @import gridExtra
#' @import stringr
#' @import viridis
#' @import RColorBrewer
#' @import BiocParallel
#' @import clusterProfiler
#' @import AnnotationHub
#' @import org.Hs.eg.db
#' @import org.Mmu.eg.db
#' @import biomaRt
#' @import SDAtools
#' @import profvis

#' 
launchHISTA <- function(...) {
  ## runApp() does not work w shiny-server
  shinyAppDir(appDir = system.file("app", package = "HISTA"))
  
}
launchHISTA.profile <- function(interval = 0.5, ...) {
  ## runApp() does not work w shiny-server
  require(profvis)
  profvis({
    shinyAppDir(appDir = system.file("app", package = "HISTA"), options = list(display.mode = "normal"))
  }, interval = interval)
}