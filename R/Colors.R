source("Help_Functions.R")

mint <- c(rgb2(228, 241, 225), rgb2(180, 217, 204), rgb2(137, 192, 182), rgb2(99, 166, 160), rgb2(68, 140, 138), rgb2(40, 114, 116), rgb2
          (13, 88, 95))
blugrn <- c(rgb2(196, 230, 195), rgb2(150, 210, 164), rgb2(109, 188, 144), rgb2(77, 162, 132), rgb2(54, 135, 122), rgb2(38, 107, 110), rgb2(29, 79, 96))
darkmint <- c(rgb2(210, 251, 212), rgb2(165, 219, 194), rgb2(123, 188, 176), rgb2(85, 156, 158), rgb2(58, 124, 137), rgb2(35, 93, 114), rgb2(18, 63, 90))

ggblue <- scales::seq_gradient_pal("#132B43", "#56B1F7", "Lab")(seq(0,1,length.out=10))

ownColors1 <- cs1 <- c("#fff7fb","#ece7f2","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#0570b0","#045a8d","#023858")

adobeColorsDiscrete <- c("#323E40","#F2AB27","#BF6B04","#732002","#D95323","#254021","#2E5902","#024873")

reverseColorSheme <- function(colors){
  n <- length(colors)
  return(colors[n:1])
}