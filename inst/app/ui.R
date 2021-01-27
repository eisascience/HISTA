# .libPaths(c("/home/groups/monkeydo/R_LIBS/3.6.1_em"))
# setwd("/home/groups/monkeydo/acc/shiny/apps/HISTA")

library(shinydashboard)

dashboardPage(
  dashboardHeader(),
  dashboardHeader(dropdownMenuOutput("messageMenu")),
  dashboardSidebar(),
  dashboardBody()
)
