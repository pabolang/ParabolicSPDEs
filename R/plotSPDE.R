#' Quick-plot SPDE sample data
#'
#' Plotting a SPDE solution field on a discrete grid using [plotly::plot_ly].
#' @param spdeData a numeric \code{NxM} matrix, where \code{N} denotes the temporal and \code{M} the spatial resolution of the equidistant grid.
#' @keywords Quick-plot samples of a SPDE model.
#'
#' @export
#' @return a plotly graphic of the SPDE model.
#' @seealso [ParabolicSPDEs::simulateSPDEmodel].
#' @examples
#' theta0 = 0
#' theta1 = 1
#' theta2 = 1
#' sigma = 0.5
#' numberSpatialPoints = 10
#' numberTemporalPoints = 1000
#' spde <- simulateSPDEmodel(theta0,theta1,theta2,sigma,numberSpatialPoints,numberTemporalPoints)
#'
#' plotSPDE(spde)






plotSPDE <- function(spdeData){
  require(plotly)
  numberSpatialPoints <- dim(spdeData)[2]-1
  numberTemporalPoints <- dim(spdeData)[1]-1
  y <- seq(0,1,1/numberSpatialPoints)
  t <- seq(0,1,1/numberTemporalPoints)
  plot_ly(x = ~y, y = ~t, z = ~spdeData,width = 1000, height = 1000) %>%
    add_surface(type="scatter3d",
                colorscale = "Cividis",
                cauto = F,
                cmin = min(spdeData),
                cmax = max(spdeData),
                contours = list(
                  z = list(
                    show=TRUE,
                    usecolormap=TRUE,
                    highlightcolor="#ff0000",
                    project=list(z=TRUE)
                  )
                ),
                showscale=FALSE,
                lighting = list(diffuse = 3),
                opacity = 1
    ) %>%
    layout(scene = list(aspectmode = "manual", aspectratio = list(x=1, y=1, z=1),
                        xaxis = list(title = list(text ='Space')),
                        yaxis = list(title = list(text ='Time')),
                        zaxis = list(title = list(text =''))
                        )
           )

}

