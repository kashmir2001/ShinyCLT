library(shiny)
library(ContourFunctions) #for multicolored titles
# install.packages("lars")
library(lars)

ui <- fluidPage(
  
  # App title ----
  titlePanel("DECKS: A Deeper Exploration of the Central Limit Theorem through the K-S Statistic"),
  
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output
      tabsetPanel(type = "tabs",
                   tabPanel("The Basics",
                            h3("Background"),
                            HTML(paste(
                              "The normal model, also known as the Gaussian model or the bell curve, is fundamental to 
                              statistics since it approximates the distribution of various phenomena across many disciplines. 
                              At the heart of this wonder is the infamous Central Limit Theorem (CLT) of which there are many 
                              versions. However, in classrooms and textbooks, not much explanation is given besides the statement 
                              of the CLT, perhaps along with its proof. This app aims to provide a greater understanding of why 
                              the normal model is indeed so normal in our world by simulating the inner workings of the classic CLT
                              (which is what you are probably familiar with).")),
                            h3("Some Important Properties of the Normal Model"),
                            HTML(paste("1. The normal model is unimodal, symmetric, and centered at the mean with thin tails that continue forever",
                                       "<br>2. According to the Empirical Rule, 68% of the data in a normal distribution lies within one standard deviation",
                                       "from the mean, 95% of the data lies within two standard deviations from the mean, and 99.7% of the data lies within",
                                       "three standard deviations from the mean")),
                            h3("Probability Distributions and the K-S Statistic"),
                            HTML(paste("This app uses a range of probability distributions, functions that give probabilities of possible values for a random variable.",
                                       "<br>In these simulations, the Kolmogorov-Smirnov statistic is used to measure the greatest difference between a 
                                        simulated distribution’s cumulative distribution function (also known as CDFs, they give the probability that a random variable with 
                                       a given probability distribution will be found at a value less than or equal to x) and a normal distribution’s CDF 
                                       with the mean and variance of the simulated distribution (with the exception of the Cauchy distribution, which has an undefined mean and variance-
                                       instead we compared the Cauchy distribution's CDF with a CDF of a normal distribution with the same location and scale). Essentially, a higher K-S statistic here corresponds to a 
                                       distribution that deviates more from normality. Additionally, the K-S statistic ranges from 0 to 1, where 0
                                       means the simulated CDF is the same as the normal CDF and 1 means the CDFs are entirely distinct.",
                                       "<br>Explore below by drawing from the distributions used in the app! Note that the blue curve is the theoretical 
                                       normal curve, and in the CDF plot on the left, 500 random draws are taken in the simulated black curves,
                                       while in the CDF plot on the right, the black curve is the theoretical CDF of the distribution you're exploring.
                                       These demonstrate the sample K-S statistic (since random draws are taken) compared to the
                                       theoretical, or expected, K-S statistic.")),
                    
                            br(), br(),
                            actionLink("toggleKS", "See more on the K-S statistic"),
                            conditionalPanel(
                              condition = "input.toggleKS % 2 == 1",
                              HTML(paste(
                                "In these simulations, the <a href='https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test' target='_blank'>K-S statistic</a> is used to measure the greatest difference between a 
                                 simulated distribution’s CDF and a normal distribution’s CDF 
                                 with the same mean and variance, but it has much wider applications; click the hyperlinked Wikipedia page to learn more."
                              ))
                            ),
                            #distribution toggles and graphs
                            h4("Normal distribution"), 
                            textInput("normm","Enter mean","0"),
                            textInput("norms","Enter standard deviation","1"),
                            plotOutput("normal"), #histogram
                            fluidRow(
                              splitLayout(cellWidths=c("50%","50%"),plotOutput("normalc"),plotOutput("normalt")) #CDF plots
                              ),
                            h4("Negative binomial distribution"),
                            textInput("nbinomr","Enter number of successes (>0)","3"),
                            textInput("nbinomp","Enter probability of success [0,1]",".5"),
                            plotOutput("nbinom"),
                            fluidRow(
                              splitLayout(cellWidths=c("50%","50%"), plotOutput("nbinomc"),plotOutput("nbinomt"))
                            ),
                            h4("Beta distribution"),
                            textInput("betaa","Enter α (>0)",".5"),
                            textInput("betab","Enter β (>0)",".5"),
                            plotOutput("beta"),
                            fluidRow(
                              splitLayout(cellWidths=c("50%","50%"), plotOutput("betac"),plotOutput("betat"))
                            ),
                            h4("Uniform distribution"),
                            textInput("unia","Enter a","0"),
                            textInput("unib","Enter b (>a)","1"),
                            plotOutput("uniform"),
                            fluidRow(
                              splitLayout(cellWidths=c("50%","50%"), plotOutput("uniformc"),plotOutput("uniformt"))
                            ),
                            h4("Poisson distribution"),
                            textInput("lambda","Enter λ (>0)","1"),
                            plotOutput("pois"),
                            fluidRow(
                              splitLayout(cellWidths=c("50%","50%"), plotOutput("poisc"),plotOutput("poist"))
                            ),
                            h4("Geometric distribution"),
                            textInput("geomp","Enter probability of success [0,1]",".5"),
                            plotOutput("geom"),
                            fluidRow(
                              splitLayout(cellWidths=c("50%","50%"), plotOutput("geomc"),plotOutput("geomt"))
                            ),
                            h4("Exponential distribution"),
                            textInput("lambda1","Enter λ (>0)","1"),
                            plotOutput("exp"),
                            fluidRow(
                              splitLayout(cellWidths=c("50%","50%"), plotOutput("expc"),plotOutput("expt"))
                            ),
                            h4("Hypergeometric distribution"),
                            textInput("hyperw","Enter number of white balls (w) in urn","1"),
                            textInput("hyperb","Enter number of black balls (b) in urn","1"),
                            textInput("hypertotal","Enter number of balls drawn from urn (w+b ≥ integer ≥ 0)","1"),
                            plotOutput("hyper"),
                            fluidRow(
                              splitLayout(cellWidths=c("50%","50%"), plotOutput("hyperc"),plotOutput("hypert"))
                            ),
                            h4("Binomial distribution"),
                            textInput("binn","Enter number of trials (integer ≥ 0)","10"),
                            textInput("binp","Enter probability of success [0,1]",".5"),
                            plotOutput("bin"),
                            fluidRow(
                              splitLayout(cellWidths=c("50%","50%"), plotOutput("binc"),plotOutput("bint"))
                            ),
                            h4("Cauchy distribution"),
                            textInput("loc","Enter location","0"),
                            textInput("scale","Enter scale","1"),
                            plotOutput("cauchy"),
                            fluidRow(
                              splitLayout(cellWidths=c("50%","50%"), plotOutput("cauchyc"),plotOutput("cauchyt"))
                            )),
                  tabPanel("Classic CLT", 
                           h3("Lindeberg-Levy or Classic CLT"),
                           HTML(paste("The sum and mean of many identical independently distributed (iid) 
                                      variables converge to approximately normal distributions no matter 
                                      how the variables themselves are distributed (as long as their mean 
                                      and variance are finite)")),
                           h3("Classic CLT Simulation"),
                           HTML(paste("The simulation below demonstrates the Classic CLT by showing the means of iid variables
                                      converging to approximately normal distributions. However, one distribution does not follow this pattern because
                                      it violates an assumption of the CLT- can you identify the distribution and the assumption?")),
                           div(style = "height:50px"), #gives space between text and the panel below
                           sidebarPanel(
                            textInput("seed","Enter seed","1"),
                            selectInput("simd","Pick a distribution",choice=c("Normal","Negative binomial","Beta",
                                                                               "Uniform","Poisson","Geometric",
                                                                               "Exponential","Hypergeometric",
                                                                               "Binomial","Cauchy") ),
                            uiOutput("select10"), #1st parameter
                            uiOutput("select20"), #2nd parameter
                            uiOutput("select30"), #3rd parameter
                            textInput("n","Enter sample size","10"),
                            textInput("N","Enter number of samples (≥3)","10"), #greater than 3 since 3 samples are automatically displayed which wouldn't make sense if they chose less than 3 samples
                            selectInput("axis", "Histogram of means has same x-axis as sample histograms?", choice=c("Yes","No"), selected = "Yes")),
                  
                           fluidRow(
                             splitLayout(
                                           plotOutput("classicclt"),    #plots of 3 samples
                                           plotOutput("classicclt2"),
                                           plotOutput("classicclt3"))
                             ),
                           
                           fluidRow(
                             column(align="center", width=12, offset=4,
                                    HTML(paste("The plots above show 3 samples of your chosen distribution and parameters. The vertical line in each plot highlights the mean of each sample."))
                             )
                           ),
                          
                           fluidRow(
                             column(align="center",
                                    plotOutput("classiccltmeans"),width=12,offset=4. #centers plot of means underneath plots of samples
                                    )
                             ),
                           
                           fluidRow(
                             column(align="center", width=12, offset=4,
                                    HTML(paste("The plot above shows a histogram of the means of the samples of your chosen distribution and parameters. When you adjust the x-axis of this plot to match the x-axis of the 3 sample plots above,
                                       you can easily compare the spread of the sample and the means- what do you notice? Also observe that each of the means in the 3 sample plots are highlighted here
                                      in the same color and an approximated normal curve in black is drawn- see how closely the histogram fits the normal curve when you
                                      increase the sample size."))
                             )
                           ),
                           
                           h3("Classic CLT K-S Statistic Simulation"),
                           HTML(paste("The simulation below computes the theoretical K-S statistic of your chosen distribution and parameters,
                                      as well as simulated K-S statistics based on 10,000 samples with various sample sizes. Note that the x-axis
                                      of the plot denotes the natural log of the sample sizes (as this is the default logarithm used in statistics)
                                      and the table resets after 7 columns.")),
                           div(style = "height:50px"), #gives space between text and the panel below
                           sidebarPanel(
                             selectInput("cltd","Pick a distribution",choice=c("Normal","Negative binomial","Beta",
                                                                               "Uniform","Poisson","Geometric",
                                                                               "Exponential","Hypergeometric",
                                                                               "Binomial","Cauchy") ),
                            uiOutput("select1"), #1st parameter
                            uiOutput("select2"), #1st parameter
                            uiOutput("select3"), #1st parameter
                            textInput("y","Choose max value on y-axis","1"),
                            actionButton("update", "Update table and graph")),
                           h4("Theoretical K-S Statistic and Simulated K-S Statistics for Multiple Sample Sizes"),
                           tableOutput("classict"), #table of K-S statistics
                           fluidRow(
                             column(align="center",
                                    plotOutput("classicplot"),width=12,offset=4 #centers plot of K-S statistics under table
                             )
                           ),
                           div(style = "height:250px"), #give space between plot and note
                           h5("Note on independence"),
                           HTML(paste("The CLT does not generally apply to sampling from a finite population, since there is no independence unless",
                                     "replacement is executed. However, the Finite Population Correction (FPC) solves this issue: if the",
                                     "sample size (n) is greater than 5% of the population of interest (N) then the standard error must have an additional factor",
                                     "&radic;<sup>(N - n)</sup>/<sub>(N - 1)</sub>",
                                     "since the means of the samples will not vary much if the samples are a large proportion of the",
                                     "population. On the other hand, if N is much greater than n, observe that the FPC is essentially 1, so has no effect",
                                     "on the standard error", collapse = " "))) #the collapse argument ensures there aren't line breaks between the strings
                 # tabPanel("Lyapunov CLT")
      )
    )
)
server<-function(input, output){
  #here we have for each distribution its histogram then its 2 CDF plots - reference the comments in the normal distribution histogram and CDFs 
  #for explanation on the other distributions' histograms and CDF plots
  
  #histogram
  output$normal<-renderPlot({
    set.seed(1) #seed so it matches its CDF plots
    hist(rnorm(500,as.numeric(input$normm),as.numeric(input$norms)),main="Simulated normal distribution",xlab="",ylab="Probability",col="cornflowerblue",prob=T)
  })
  
  #CDF plot with simulated CDF of chosen distribution
  output$normalc<-renderPlot({
    set.seed(1) 
    n<-500
    m<-as.numeric(input$normm)
    s<-as.numeric(input$norms)
    norm<-rnorm(n,m,s) #simulated distribution 
    
    ks_result <- ks.test(norm, "pnorm", m, s) 
    ks <- ks_result$statistic 
    ks<-round(ks,3) #KS statistic

    myecdf <- ecdf(norm) #simulated distribution's CDF
    seq <- seq(min(norm), max(norm), length.out = n) #x-axis for the normal CDF we're comparing our simulation against 
    
    plot(myecdf,xlab="",ylab="Probability",main="") #plot simulated distribution's CDF 
    lines(seq,pnorm(seq,m,s),col="deepskyblue3",lwd=5) #plot theoretical normal CDF
    
    multicolor.title(c("Simulated Normal CDF and ","Normal CDF ", "with", " K-S ",ks), c("black","deepskyblue3","black","red","red")) #colors correspond to features on graph
    
    abs_diff <- abs(myecdf(seq) - pnorm(seq,m, s)) #create vector of absolute value differences between the CDFs
    max <- which.max(abs_diff) #find max difference - that's the KS statistic
    y0 <- myecdf(seq)[max[1]] #y-value of simulated CDF where there is the greatest diff between the CDFs
    y1 <- pnorm(seq,m, s)[max[1]] #y-value of theoretical CDF where there is the greatest diff between the CDFs
    segments(seq[max[1]],y0,seq[max[1]],y1,col="red",lwd=5) #create red line to represent the KS statistic
  })
  
  #CDF plot with theoretical CDF of chosen distribution
  output$normalt<-renderPlot({
    set.seed(1)
    n<-500
    m<-as.numeric(input$normm)
    s<-as.numeric(input$norms)
    norm<-rnorm(n,m,s)
  
    seq <- seq(min(norm), max(norm), length.out = n)
    ncdf<-pnorm(seq,m,s) #theoretical CDF of chosen distribution 
    normcdf<-pnorm(seq,m,s) #theoretical normal CDF we're comparing against 
    
    diff<-abs(ncdf-normcdf) #create vector of absolute value differences between the CDFs
    ks<-max(diff) #theoretical K-S stat
    
    plot(seq,pnorm(seq,m,s),xlab="",ylab="Probability",main="",pch=1,cex=1,ylim=c(0,1)) #plot simulated distribution's theoretical CDF 
    lines(seq,pnorm(seq,m,s),col="deepskyblue3",lwd=5) #plot the same theoretical normal distribution's CDF we're comparing against
    
    multicolor.title(c("Normal CDF and ","Normal CDF ", "with", " K-S ",ks), c("black","deepskyblue3","black","red","red"))
    
    abs_diff <- abs(pnorm(seq,m,s) - pnorm(seq,m, s)) #create vector of absolute value differences between the CDFs
    max <- which.max(abs_diff) #find max difference - that's the KS statistic
    y0 <- pnorm(seq,m,s)[max[1]] #y-value of theoretical CDF of chosen distribution where there is the greatest diff between the CDFs
    y1 <- pnorm(seq,m, s)[max[1]] #y-value of the normal theoretical CDF we're comparing against where there is the greatest diff between the CDFs
    segments(seq[max[1]],y0,seq[max[1]],y1,col="red",lwd=5) #create red line to represent the KS statistic
  })
  output$nbinom<-renderPlot({
    set.seed(1)
    r<-as.numeric(input$nbinomr)
    p<-as.numeric(input$nbinomp)
    n<-500
    nbinom<-rnbinom(n,r,p)
    hist(nbinom,main="Simulated negative binomial distribution",xlab="",ylab="Probability",col="cornflowerblue",prob=T)
  })
  output$nbinomc<-renderPlot({
    set.seed(1)
    n<-500
    r<-as.numeric(input$nbinomr)
    p<-as.numeric(input$nbinomp)
    nbinom<-rnbinom(n,r,p)
    
    ks_result <- ks.test(nbinom, "pnorm", (r*(1-p))/p, sqrt((r*(1-p))/(p^2))) #use mean and sd of negative binomial distribution as mean/sd of normal distribution we're comparing against
    ks <- ks_result$statistic
    ks<-round(ks,3)
    
    myecdf <- ecdf(nbinom)
    seq <- seq(min(nbinom), max(nbinom), length.out = n)

    plot(myecdf,xlab="",ylab="Probability",main="")
    lines(seq,pnorm(seq,(r*(1-p))/p, sqrt((r*(1-p))/(p^2))),col="deepskyblue3",lwd=5)
    
    multicolor.title(c("Simulated Negative Binomial CDF and ","Normal CDF ", "with", " K-S ",ks), c("black","deepskyblue3","black","red","red"))
    
    abs_diff <- abs(myecdf(seq) - pnorm(seq,(r*(1-p))/p, sqrt((r*(1-p))/(p^2)))) 
    max <- which.max(abs_diff)
    y0 <- myecdf(seq)[max[1]]
    y1 <- pnorm(seq,(r*(1-p))/p, sqrt((r*(1-p))/(p^2)))[max[1]]
    segments(seq[max[1]],y0,seq[max[1]],y1,col="red",lwd=5)
  })
  output$nbinomt<-renderPlot({
    set.seed(1)
    n<-500
    r<-as.numeric(input$nbinomr)
    p<-as.numeric(input$nbinomp)
    nbinom<-rnbinom(n,r,p)
    
    seq <- seq(min(nbinom), max(nbinom), length.out = n)
    cdf<-pnbinom(seq,r,p) #theoretical cdf of chosen distribution
    normcdf<-pnorm(seq,(r*(1-p))/p,sqrt((r*(1-p))/(p^2))) #normal CDF
    
    diff<-abs(cdf-normcdf)
    ks<-max(diff) #theoretical K-S stat
    ks<-round(ks,3)
    
    plot(seq,cdf,xlab="",ylab="Probability",main="",pch=1,cex=1,ylim=c(0,1))
    lines(seq,normcdf,col="deepskyblue3",lwd=5)
    
    multicolor.title(c("Negative Binomial CDF and ","Normal CDF ", "with", " K-S ",ks), c("black","deepskyblue3","black","red","red"))
    
    abs_diff <- abs(cdf - normcdf) 
    max <- which.max(abs_diff)
    y0 <- cdf[max[1]]
    y1 <- normcdf[max[1]]
    segments(seq[max[1]],y0,seq[max[1]],y1,col="red",lwd=5)
  })
  output$beta<-renderPlot({
    set.seed(1)
    a<-as.numeric(input$betaa)
    b<-as.numeric(input$betab)
    hist(rbeta(500,a,b),main="Simulated beta distribution",xlab="",ylab="Probability",col="cornflowerblue",prob=T)
  })
  output$betac<-renderPlot({
    set.seed(1)
    n<-500
    a<-as.numeric(input$betaa)
    b<-as.numeric(input$betab)
    beta<-rbeta(n,a,b)
    
    ks_result <- ks.test(beta, "pnorm", a/(a+b), sqrt((a*b)/((a+b)^2*(a+b+1))))
    ks <- ks_result$statistic
    ks<-round(ks,3)
    
    myecdf <- ecdf(beta)
    seq <- seq(min(beta), max(beta), length.out = n)
    
    plot(myecdf,xlab="",ylab="Probability",main="")
    lines(seq,pnorm(seq,a/(a+b), sqrt((a*b)/((a+b)^2*(a+b+1)))),col="deepskyblue3",lwd=5)
    
    multicolor.title(c("Simulated Beta CDF and ","Normal CDF ", "with", " K-S ",ks), c("black","deepskyblue3","black","red","red"))
    
    abs_diff <- abs(myecdf(seq) - pnorm(seq,a/(a+b), sqrt((a*b)/((a+b)^2*(a+b+1))))) 
    max <- which.max(abs_diff)
    y0 <- myecdf(seq)[max[1]]
    y1 <- pnorm(seq,a/(a+b), sqrt((a*b)/((a+b)^2*(a+b+1))))[max[1]]
    segments(seq[max[1]],y0,seq[max[1]],y1,col="red",lwd=5)
  })
  output$betat<-renderPlot({
    set.seed(1)
    n<-500
    a<-as.numeric(input$betaa)
    b<-as.numeric(input$betab)
    beta<-rbeta(n,a,b)
    
    seq <- seq(min(beta), max(beta), length.out = n)
    cdf<-pbeta(seq,a,b)
    normcdf<-pnorm(seq,a/(a+b),sqrt((a*b)/((a+b)^2*(a+b+1))))
    
    diff<-abs(cdf-normcdf)
    ks<-max(diff) #theoretical K-S stat
    ks<-round(ks,3)
    
    plot(seq,cdf,xlab="",ylab="Probability",main="",pch=1,cex=1,ylim=c(0,1))
    lines(seq,normcdf,col="deepskyblue3",lwd=5)
    
    multicolor.title(c("Beta CDF and ","Normal CDF ", "with", " K-S ",ks), c("black","deepskyblue3","black","red","red"))
    
    abs_diff <- abs(cdf - normcdf) 
    max <- which.max(abs_diff)
    y0 <- cdf[max[1]]
    y1 <- normcdf[max[1]]
    segments(seq[max[1]],y0,seq[max[1]],y1,col="red",lwd=5)
  })
  output$uniform<-renderPlot({
    set.seed(1)
    hist(runif(500,as.numeric(input$unia),as.numeric(input$unib)),main="Simulated uniform distribution",xlab="",ylab="Probability",col="cornflowerblue",prob=T)
  })
  output$uniformc<-renderPlot({
    set.seed(1)
    n<-500
    a<-as.numeric(input$unia)
    b<-as.numeric(input$unib)
    uni<-runif(n,a,b)
    
    ks_result <- ks.test(uni, "pnorm", .5*(a+b), sqrt((1/12)*(b-a)^2))
    ks <- ks_result$statistic
    ks<-round(ks,3)
    
    myecdf <- ecdf(uni)
    seq <- seq(min(uni), max(uni), length.out = n)
    
    plot(myecdf,xlab="",ylab="Probability",main="")
    lines(seq,pnorm(seq,.5*(a+b), sqrt((1/12)*(b-a)^2)),col="deepskyblue3",lwd=5)
    
    multicolor.title(c("Simulated Uniform CDF and ","Normal CDF ", "with", " K-S ",ks), c("black","deepskyblue3","black","red","red"))
    
    abs_diff <- abs(myecdf(seq) - pnorm(seq,.5*(a+b), sqrt((1/12)*(b-a)^2)))
    max <- which.max(abs_diff)
    y0 <- myecdf(seq)[max[1]]
    y1 <- pnorm(seq,.5*(a+b), sqrt((1/12)*(b-a)^2))[max[1]]
    segments(seq[max[1]],y0,seq[max[1]],y1,col="red",lwd=5)
  })
  output$uniformt<-renderPlot({
    set.seed(1)
    n<-500
    a<-as.numeric(input$unia)
    b<-as.numeric(input$unib)
    uni<-runif(n,a,b)
    
    seq <- seq(min(uni), max(uni), length.out = n)
    cdf<-punif(seq,a,b)
    normcdf<-pnorm(seq,.5*(a+b),sqrt((1/12)*(b-a)^2))
    
    diff<-abs(cdf-normcdf)
    ks<-max(diff) #theoretical K-S stat
    ks<-round(ks,3)
    
    plot(seq,cdf,xlab="",ylab="Probability",main="",pch=1,cex=1,ylim=c(0,1))
    lines(seq,normcdf,col="deepskyblue3",lwd=5)
    
    multicolor.title(c("Uniform CDF and ","Normal CDF ", "with", " K-S ",ks), c("black","deepskyblue3","black","red","red"))
    
    abs_diff <- abs(cdf - normcdf) 
    max <- which.max(abs_diff)
    y0 <- cdf[max[1]]
    y1 <- normcdf[max[1]]
    segments(seq[max[1]],y0,seq[max[1]],y1,col="red",lwd=5)
  })
  output$pois<-renderPlot({
    set.seed(1)
    hist(rpois(500,as.numeric(input$lambda)),main="Simulated poisson distribution",xlab="",ylab="Probability",col="cornflowerblue",prob=T)
  })
  output$poisc<-renderPlot({
    set.seed(1)
    n<-500
    l<-as.numeric(input$lambda)
    pois<-rpois(n,l)
    
    ks_result <- ks.test(pois, "pnorm", l, sqrt(l))
    ks <- ks_result$statistic
    ks<-round(ks,3)
    
    myecdf <- ecdf(pois)
    seq <- seq(min(pois), max(pois), length.out = n)
    
    plot(myecdf,xlab="",ylab="Probability",main="")
    lines(seq,pnorm(seq,l, sqrt(l)),col="deepskyblue3",lwd=5)
    
    multicolor.title(c("Simulated Poisson CDF and ","Normal CDF ", "with", " K-S ",ks), c("black","deepskyblue3","black","red","red"))
    
    abs_diff <- abs(myecdf(seq) - pnorm(seq,l, sqrt(l))) 
    max <- which.max(abs_diff)
    y0 <- myecdf(seq)[max[1]]
    y1 <- pnorm(seq,l, sqrt(l))[max[1]]
    segments(seq[max[1]],y0,seq[max[1]],y1,col="red",lwd=5)
  })
  output$poist<-renderPlot({
    set.seed(1)
    n<-500
    l<-as.numeric(input$lambda)
    pois<-rpois(n,l)
    
    seq <- seq(min(pois), max(pois), length.out = n)
    cdf<-ppois(seq,l)
    normcdf<-pnorm(seq,l,sqrt(l))
    
    diff<-abs(cdf-normcdf)
    ks<-max(diff) #theoretical K-S stat
    ks<-round(ks,3)
    
    plot(seq,cdf,xlab="",ylab="Probability",main="",pch=1,cex=1,ylim=c(0,1))
    lines(seq,normcdf,col="deepskyblue3",lwd=5)
    
    multicolor.title(c("Poisson CDF and ","Normal CDF ", "with", " K-S ",ks), c("black","deepskyblue3","black","red","red"))
    
    abs_diff <- abs(cdf - normcdf) 
    max <- which.max(abs_diff)
    y0 <- cdf[max[1]]
    y1 <- normcdf[max[1]]
    segments(seq[max[1]],y0,seq[max[1]],y1,col="red",lwd=5)
  })
  output$geom<-renderPlot({
    set.seed(1)
    hist(rgeom(500,as.numeric(input$geomp)),main="Simulated geometric distribution",xlab="",ylab="Probability",col="cornflowerblue",prob=T)
  })
  output$geomc<-renderPlot({
    set.seed(1)
    n<-500
    p<-as.numeric(input$geomp)
    geom<-rgeom(n,p)
    
    ks_result <- ks.test(geom, "pnorm", (1-p)/p, sqrt((1-p)/p^2))
    ks <- ks_result$statistic
    ks<-round(ks,3)
    
    myecdf <- ecdf(geom)
    seq <- seq(-10, max(geom), length.out = n) #the geometric distribution is only positive, but we want to show negative values too since the normal distribution covers all real #s
    
    plot(myecdf,xlab="",ylab="Probability",main="",xlim=c(seq[1],seq[length(seq)])) #x-limits defined explicitly here to match seq so it will include negative values
    lines(seq,pnorm(seq,(1-p)/p, sqrt((1-p)/p^2)),col="deepskyblue3",lwd=5)
    
    multicolor.title(c("Simulated Geometric CDF and ","Normal CDF ", "with", " K-S ",ks), c("black","deepskyblue3","black","red","red"))
    
    abs_diff <- abs(myecdf(seq) - pnorm(seq,(1-p)/p, sqrt((1-p)/p^2))) 
    max <- which.max(abs_diff)
    y0 <- myecdf(seq)[max[1]]
    y1 <- pnorm(seq,(1-p)/p, sqrt((1-p)/p^2))[max[1]]
    segments(seq[max[1]],y0,seq[max[1]],y1,col="red",lwd=5)
  })
  output$geomt<-renderPlot({
    set.seed(1)
    n<-500
    p<-as.numeric(input$geomp)
    geom<-rgeom(n,p)
    
    seq <- seq(-10, max(geom), length.out = n) 
    cdf<-pgeom(seq,p)
    normcdf<-pnorm(seq,(1-p)/p, sqrt((1-p)/p^2))
    
    diff<-abs(cdf-normcdf)
    ks<-max(diff) #theoretical K-S stat
    ks<-round(ks,3)
    
    plot(seq,cdf,xlab="",ylab="Probability",main="",pch=1,cex=1,ylim=c(0,1))
    lines(seq,normcdf,col="deepskyblue3",lwd=5)
    
    multicolor.title(c("Geometric CDF and ","Normal CDF ", "with", " K-S ",ks), c("black","deepskyblue3","black","red","red"))
    
    abs_diff <- abs(cdf - normcdf) 
    max <- which.max(abs_diff)
    y0 <- cdf[max[1]]
    y1 <- normcdf[max[1]]
    segments(seq[max[1]],y0,seq[max[1]],y1,col="red",lwd=5)
  })
  output$exp<-renderPlot({
    set.seed(1)
    hist(rexp(500,as.numeric(input$lambda1)),main="Simulated exponential distribution",xlab="",ylab="Probability",col="cornflowerblue",prob=T)
  })
  output$expc<-renderPlot({
    set.seed(1)
    n<-500
    l<-as.numeric(input$lambda1)
    exp<-rexp(n,l)
    
    ks_result <- ks.test(exp, "pnorm", 1/l, 1/l)
    ks <- ks_result$statistic
    ks<-round(ks,3)
    
    myecdf <- ecdf(exp)
    seq <- seq(-10, max(exp), length.out = n) #the exponential distribution is only positive, but we want to show negative values too since the normal distribution covers all real #s
    
    plot(myecdf,xlab="",ylab="Probability",main="",xlim=c(seq[1],seq[length(seq)]))
    lines(seq,pnorm(seq,1/l, 1/l),col="deepskyblue3",lwd=5)
    
    multicolor.title(c("Simulated Exponential CDF and ","Normal CDF ", "with", " K-S ",ks), c("black","deepskyblue3","black","red","red"))
    
    abs_diff <- abs(myecdf(seq) - pnorm(seq,1/l, 1/l)) 
    max <- which.max(abs_diff)
    y0 <- myecdf(seq)[max[1]]
    y1 <- pnorm(seq,1/l, 1/l)[max[1]]
    segments(seq[max[1]],y0,seq[max[1]],y1,col="red",lwd=5)
  })
  output$expt<-renderPlot({
    set.seed(1)
    n<-500
    l<-as.numeric(input$lambda1)
    exp<-rexp(n,l)
    
    seq <- seq(-10, max(exp), length.out = n)
    cdf<-pexp(seq,l)
    normcdf<-pnorm(seq,1/l, 1/l)
    
    diff<-abs(cdf-normcdf)
    ks<-max(diff) #theoretical K-S stat
    ks<-round(ks,3)
    
    plot(seq,cdf,xlab="",ylab="Probability",main="",pch=1,cex=1,ylim=c(0,1))
    lines(seq,normcdf,col="deepskyblue3",lwd=5)
    
    multicolor.title(c("Exponential CDF and ","Normal CDF ", "with", " K-S ",ks), c("black","deepskyblue3","black","red","red"))
    
    abs_diff <- abs(cdf - normcdf) 
    max <- which.max(abs_diff)
    y0 <- cdf[max[1]]
    y1 <- normcdf[max[1]]
    segments(seq[max[1]],y0,seq[max[1]],y1,col="red",lwd=5)
  })
  output$hyper<-renderPlot({
    set.seed(1)
    hist(rhyper(500,as.numeric(input$hyperw),as.numeric(input$hyperb),as.numeric(input$hypertotal)),main="Simulated hypergeometric distribution",xlab="",ylab="Probability",col="cornflowerblue",prob=T)
  })
  output$hyperc<-renderPlot({
    set.seed(1)
    n<-500
    w<-as.numeric(input$hyperw)
    b<-as.numeric(input$hyperb)
    t<-as.numeric(input$hypertotal)
    hyp<-rhyper(n,w,b,t)
    
    ks_result <- ks.test(hyp, "pnorm", t*(w/(w+b)), sqrt(t*(w/(w+b))*(1-(w/(w+b)))*((w+b-t)/(w+b-1))))
    ks <- ks_result$statistic
    ks<-round(ks,3)
    
    myecdf <- ecdf(hyp)
    seq <- seq(min(hyp), max(hyp), length.out = n)
    
    plot(myecdf,xlab="",ylab="Probability",main="")
    lines(seq,pnorm(seq,t*(w/(w+b)), sqrt(t*(w/(w+b))*(1-(w/(w+b)))*((w+b-t)/(w+b-1)))),col="deepskyblue3",lwd=5)
    
    multicolor.title(c("Simulated Hypergeometric CDF and ","Normal CDF ", "with", " K-S ",ks), c("black","deepskyblue3","black","red","red"))
    
    abs_diff <- abs(myecdf(seq) - pnorm(seq,t*(w/(w+b)), sqrt(t*(w/(w+b))*(1-(w/(w+b)))*((w+b-t)/(w+b-1))))) 
    max <- which.max(abs_diff)
    y0 <- myecdf(seq)[max[1]]
    y1 <- pnorm(seq,t*(w/(w+b)), sqrt(t*(w/(w+b))*(1-(w/(w+b)))*((w+b-t)/(w+b-1))))[max[1]]
    segments(seq[max[1]],y0,seq[max[1]],y1,col="red",lwd=5)
  })
  output$hypert<-renderPlot({
    set.seed(1)
    n<-500
    w<-as.numeric(input$hyperw)
    b<-as.numeric(input$hyperb)
    t<-as.numeric(input$hypertotal)
    hyp<-rhyper(n,w,b,t)
    
    seq <- seq(min(hyp), max(hyp), length.out = n)
    cdf<-phyper(seq,w,b,t)
    normcdf<-pnorm(seq,t*(w/(w+b)), sqrt(t*(w/(w+b))*(1-(w/(w+b)))*((w+b-t)/(w+b-1))))
    
    diff<-abs(cdf-normcdf)
    ks<-max(diff) #theoretical K-S stat
    ks<-round(ks,3)
    
    plot(seq,cdf,xlab="",ylab="Probability",main="",pch=1,cex=1,ylim=c(0,1))
    lines(seq,normcdf,col="deepskyblue3",lwd=5)
    
    multicolor.title(c("Hypergeometric CDF and ","Normal CDF ", "with", " K-S ",ks), c("black","deepskyblue3","black","red","red"))
    
    abs_diff <- abs(cdf - normcdf) 
    max <- which.max(abs_diff)
    y0 <- cdf[max[1]]
    y1 <- normcdf[max[1]]
    segments(seq[max[1]],y0,seq[max[1]],y1,col="red",lwd=5)
  })
  output$bin<-renderPlot({
    set.seed(1)
    hist(rbinom(500,as.numeric(input$binn),as.numeric(input$binp)),main="Simulated binomial distribution",xlab="",ylab="Probability",col="cornflowerblue",prob=T)
  })
  output$binc<-renderPlot({
    set.seed(1)
    n<-500
    t<-as.numeric(input$binn)
    p<-as.numeric(input$binp)
    bin<-rbinom(n,t,p)
    
    ks_result <- ks.test(bin, "pnorm", t*p, sqrt(t*p*(1-p)))
    ks <- ks_result$statistic
    ks<-round(ks,3)
    
    myecdf <- ecdf(bin)
    seq <- seq(min(bin), max(bin), length.out = n)
    
    plot(myecdf,xlab="",ylab="Probability",main="")
    lines(seq,pnorm(seq,t*p, sqrt(t*p*(1-p))),col="deepskyblue3",lwd=5)
    
    multicolor.title(c("Simulated Binomial CDF and ","Normal CDF ", "with", " K-S ",ks), c("black","deepskyblue3","black","red","red"))
    
    abs_diff <- abs(myecdf(seq) - pnorm(seq,t*p, sqrt(t*p*(1-p)))) 
    max <- which.max(abs_diff)
    y0 <- myecdf(seq)[max[1]]
    y1 <- pnorm(seq,t*p, sqrt(t*p*(1-p)))[max[1]]
    segments(seq[max[1]],y0,seq[max[1]],y1,col="red",lwd=5)
  })
  output$bint<-renderPlot({
    set.seed(1)
    n<-500
    t<-as.numeric(input$binn)
    p<-as.numeric(input$binp)
    bin<-rbinom(n,t,p)
    
    seq <- seq(min(bin), max(bin), length.out = n)
    cdf<-pbinom(seq,t,p)
    normcdf<-pnorm(seq,t*p, sqrt(t*p*(1-p)))
    
    diff<-abs(cdf-normcdf)
    ks<-max(diff) #theoretical K-S stat
    ks<-round(ks,3)
    
    plot(seq,cdf,xlab="",ylab="Probability",main="",pch=1,cex=1,ylim=c(0,1))
    lines(seq,normcdf,col="deepskyblue3",lwd=5)
    
    multicolor.title(c("Binomial CDF and ","Normal CDF ", "with", " K-S ",ks), c("black","deepskyblue3","black","red","red"))
    
    abs_diff <- abs(cdf - normcdf) 
    max <- which.max(abs_diff)
    y0 <- cdf[max[1]]
    y1 <- normcdf[max[1]]
    segments(seq[max[1]],y0,seq[max[1]],y1,col="red",lwd=5)
  })
  output$cauchy<-renderPlot({
    set.seed(1)
    hist(rcauchy(500,as.numeric(input$loc),as.numeric(input$scale)),main="Simulated cauchy distribution",xlab="",ylab="Probability",col="cornflowerblue",prob=T)
  })
  output$cauchyc<-renderPlot({
    set.seed(1)
    n<-500
    l<-as.numeric(input$loc)
    s<-as.numeric(input$scale)
    cauchy<-rcauchy(n,l,s)
    
    #In my independent study in spring '22, we estimated Cauchy's mean and variance through the forloop below, where we took the average of many 
    #Cauchy samples' means and variances in order to give the Normal distribution which we're comparing against for the K-S statistic a mean
    #and variance. However, we are now just assigning the Normal distribution the same location and scale as the Cauchy distribution, since
    #Cauchy technically doesn't have a mean/variance.
    #sims<-1000
    # storagem<-rep(NA,sims) 
    # storagev<-rep(NA,sims)
    # for (ii in 1:sims) {
    #   storagem[ii]<-mean(rcauchy(n))
    #   storagev[ii]<-var(rcauchy(n))
    #   mean<-mean(storagem)
    #   var<-mean(storagev)
    # }
    ks=ks.test(cauchy, "pnorm", l, s)$statistic
    ks<-round(ks,3)
    
    myecdf <- ecdf(cauchy)
    seq <- seq(min(cauchy), max(cauchy), length.out = n)
    
    plot(myecdf,xlab="",ylab="Probability",main="")
    lines(seq,pnorm(seq,l, s),col="deepskyblue3",lwd=5)
    
    multicolor.title(c("Simulated Cauchy CDF and ","Normal CDF ", "with", " K-S ",ks), c("black","deepskyblue3","black","red","red"))
    
    abs_diff <- abs(myecdf(seq) - pnorm(seq,l, s)) 
    max <- which.max(abs_diff)
    y0 <- myecdf(seq)[max[1]]
    y1 <- pnorm(seq,l, s)[max[1]]
    segments(seq[max[1]],y0,seq[max[1]],y1,col="red",lwd=5)
  })
  output$cauchyt<-renderPlot({
    set.seed(1)
    n<-500
    l<-as.numeric(input$loc)
    s<-as.numeric(input$scale)
    cauchy<-rcauchy(n,l,s)
    
    seq <- seq(min(cauchy), max(cauchy), length.out = n)
    cdf<-pcauchy(seq,l,s)
    normcdf<-pnorm(seq,l, s)
    
    diff<-abs(cdf-normcdf)
    ks<-max(diff) #theoretical K-S stat
    ks<-round(ks,3)
    
    plot(seq,cdf,xlab="",ylab="Probability",main="",pch=1,cex=1,ylim=c(0,1))
    lines(seq,normcdf,col="deepskyblue3",lwd=5)
    
    multicolor.title(c("Cauchy CDF and ","Normal CDF ", "with", " K-S ",ks), c("black","deepskyblue3","black","red","red"))
    
    abs_diff <- abs(cdf - normcdf) 
    max <- which.max(abs_diff)
    y0 <- cdf[max[1]]
    y1 <- normcdf[max[1]]
    segments(seq[max[1]],y0,seq[max[1]],y1,col="red",lwd=5)
  })
  output$select10 <- renderUI({                  #Based on the distribution you chose, it gives that distribution's first parameter
    if (input$simd == "Normal") {
      textInput("norm1", "Enter mean","0")
    }
    else if (input$simd == "Negative binomial") {
      textInput("nb1", "Enter number of successes (>0)","3")
    }
    else if (input$simd == "Beta") {
      textInput("beta1", "Enter α (>0)",".5")
    }
    else if (input$simd == "Uniform") {
      textInput("unif1", "Enter a","0")
    }
    else if (input$simd == "Poisson") {
      textInput("pois1", "Enter λ (>0)","1")
    }
    else if (input$simd == "Geometric") {
      textInput("geom1", "Enter probability of success [0,1]",".5")
    }
    else if (input$simd == "Exponential") {
      textInput("exp1", "Enter λ (>0)","1")
    }
    else if (input$simd == "Hypergeometric") {
      textInput("hyp1", "Enter number of white balls (w) in urn","1")
    }
    else if (input$simd == "Binomial") {
      textInput("bin1", "Enter number of trials (integer ≥ 0)","10")
    }
    else {
      textInput("cauchy1", "Enter location","0")
    }
  })
  output$select20 <- renderUI({                   #Based on the distribution you chose, it gives that distribution's second parameter         
    if (input$simd == "Normal") {
      textInput("norm2", "Enter standard deviation","1")
    }
    else if (input$simd == "Negative binomial") {
      textInput("nb2", "Enter probability of success [0,1]",".5")
    }
    else if (input$simd == "Beta") {
      textInput("beta2", "Enter β (>0)",".5")
    }
    else if (input$simd == "Uniform") {
      textInput("unif2", "Enter b (>a)","1")
    }
    else if (input$simd == "Hypergeometric") {
      textInput("hyp2", "Enter number of black balls (b) in urn","1")
    }
    else if (input$simd == "Binomial") {
      textInput("bin2", "Enter probability of success [0,1]",".5")
    }
    else if (input$simd == "Cauchy"){
      textInput("cauchy2", "Enter scale","1")
    }
    else {
    }
  })
  output$select30 <- renderUI({                      #Based on the distribution you chose, it gives that distribution's third parameter
    if (input$simd == "Hypergeometric") {
      textInput("hyp3", "Enter number of balls drawn from urn (w+b ≥ integer ≥ 0)","1")
    }
    else {
    }
  })
  
  #create objects that will store reactive values - in this case, it will store 3 samples of the distribution/parameters the user chose
  first_data <- reactiveVal()
  second_data <- reactiveVal()
  third_data <- reactiveVal()
  
  #assign those objects values based on the distribution/parameters user chose
  observe({
    if (input$simd == "Normal") {
      set.seed(input$seed)
      first <- rnorm(as.numeric(input$n), as.numeric(input$norm1), as.numeric(input$norm2))
      first_data(first)  # Store the value of 'first' in the reactiveVal
      second <- rnorm(as.numeric(input$n), as.numeric(input$norm1), as.numeric(input$norm2))
      second_data(second)  
      third <- rnorm(as.numeric(input$n), as.numeric(input$norm1), as.numeric(input$norm2))
      third_data(third)  
    }
    else if (input$simd == "Negative binomial") {
      set.seed(input$seed)
      first <- rnbinom(as.numeric(input$n), as.numeric(input$nb1), as.numeric(input$nb2))
      first_data(first)  # Store the value of 'first' in the reactiveVal
      second <- rnbinom(as.numeric(input$n), as.numeric(input$nb1), as.numeric(input$nb2))
      second_data(second)  
      third <- rnbinom(as.numeric(input$n), as.numeric(input$nb1), as.numeric(input$nb2))
      third_data(third)  
    }
    else if (input$simd == "Beta") {
      set.seed(input$seed)
      first <- rbeta(as.numeric(input$n), as.numeric(input$beta1), as.numeric(input$beta2))
      first_data(first)  # Store the value of 'first' in the reactiveVal
      second <- rbeta(as.numeric(input$n), as.numeric(input$beta1), as.numeric(input$beta2))
      second_data(second)  
      third <- rbeta(as.numeric(input$n), as.numeric(input$beta1), as.numeric(input$beta2))
      third_data(third)  
    }
    else if (input$simd == "Uniform") {
      set.seed(input$seed)
      first <- runif(as.numeric(input$n), as.numeric(input$unif1), as.numeric(input$unif2))
      first_data(first)  # Store the value of 'first' in the reactiveVal
      second <- runif(as.numeric(input$n), as.numeric(input$unif1), as.numeric(input$unif2))
      second_data(second)  
      third <- runif(as.numeric(input$n), as.numeric(input$unif1), as.numeric(input$unif2))
      third_data(third)  
    }
    else if (input$simd == "Poisson") {
      set.seed(input$seed)
      first <- rpois(as.numeric(input$n), as.numeric(input$pois1))
      first_data(first)  
      second <- rpois(as.numeric(input$n), as.numeric(input$pois1))
      second_data(second)  
      third <- rpois(as.numeric(input$n), as.numeric(input$pois1))
      third_data(third)  
    }
    else if (input$simd == "Geometric") {
      set.seed(input$seed)
      first <- rgeom(as.numeric(input$n), as.numeric(input$geom1))
      first_data(first)  
      second <- rgeom(as.numeric(input$n), as.numeric(input$geom1))
      second_data(second)  
      third <- rgeom(as.numeric(input$n), as.numeric(input$geom1))
      third_data(third)  
    }
    else if (input$simd == "Exponential") {
      set.seed(input$seed)
      first <- rexp(as.numeric(input$n), as.numeric(input$exp1))
      first_data(first)  
      second <- rexp(as.numeric(input$n), as.numeric(input$exp1))
      second_data(second)  
      third <- rexp(as.numeric(input$n), as.numeric(input$exp1))
      third_data(third)  
    }
    else if (input$simd == "Hypergeometric") {
      set.seed(input$seed)
      first <- rhyper(as.numeric(input$n), as.numeric(input$hyp1),as.numeric(input$hyp2),as.numeric(input$hyp3))
      first_data(first)  
      second <- rhyper(as.numeric(input$n), as.numeric(input$hyp1),as.numeric(input$hyp2),as.numeric(input$hyp3))
      second_data(second)  
      third <- rhyper(as.numeric(input$n), as.numeric(input$hyp1),as.numeric(input$hyp2),as.numeric(input$hyp3))
      third_data(third)  
    }
    else if (input$simd == "Binomial") {
      set.seed(input$seed)
      first <- rbinom(as.numeric(input$n), as.numeric(input$bin1),as.numeric(input$bin2))
      first_data(first)  
      second <- rbinom(as.numeric(input$n), as.numeric(input$bin1),as.numeric(input$bin2))
      second_data(second)  
      third <- rbinom(as.numeric(input$n), as.numeric(input$bin1),as.numeric(input$bin2))
      third_data(third)  
    }
    else {
      set.seed(input$seed)
      first <- rcauchy(as.numeric(input$n), as.numeric(input$cauchy1),as.numeric(input$cauchy2))
      first_data(first)  
      second <- rcauchy(as.numeric(input$n), as.numeric(input$cauchy1),as.numeric(input$cauchy2))
      second_data(second)  
      third <- rcauchy(as.numeric(input$n), as.numeric(input$cauchy1),as.numeric(input$cauchy2))
      third_data(third)  
    }
  })
  
  #plot 1st sample
  output$classicclt <- renderPlot({
    data1 <- first_data()  # Access the value of 'first' from the reactiveVal
    data2 <- second_data()  
    data3 <- third_data()  
    xmin<-min(data1,data2,data3) #we want all the samples to have the same x-axis so they're comparable, so we'll take the smallest and largest values from all samples as the limits
    xmax<-max(data1,data2,data3)
    mean<-round(mean(data1),3)
    
    hist(data1, xlab = "", col = "magenta4", main = "",xlim=c(floor(xmin),ceiling(xmax))) #we want the limits to be integers (so it's clean) that still capture the smallest and largest values so we use floor and ceiling
    abline(v=mean,col="mediumorchid1",lwd=5) #highlight mean
    multicolor.title(c("1st sample with ","sample mean ", mean),c("magenta4","mediumorchid1","mediumorchid1"))

  })
  
  #plot 2nd sample
  output$classicclt2<-renderPlot({

    data1 <- first_data()  # Access the value of 'first' from the reactiveVal
    data2 <- second_data()  
    data3 <- third_data()  
    xmin<-min(data1,data2,data3)
    xmax<-max(data1,data2,data3)
    mean<-round(mean(data2),3)

    hist(data2,xlab="",col="maroon1",main="",xlim=c(floor(xmin),ceiling(xmax)))
    abline(v=mean,col="salmon1",lwd=5)
    h3(multicolor.title(c("2nd sample with ","sample mean ", mean),c("maroon1","salmon1","salmon1")))

  })
  
  #plot 3rd sample
  output$classicclt3<-renderPlot({
    data1 <- first_data()  # Access the value of 'first' from the reactiveVal
    data2 <- second_data()  
    data3 <- third_data()  
    xmin<-min(data1,data2,data3)
    xmax<-max(data1,data2,data3)
    mean<-round(mean(data3),3)
    
    hist(data3,xlab="",col="darkgreen",main="",xlim=c(floor(xmin),ceiling(xmax)))
    abline(v=mean,col="green3",lwd=5)
    multicolor.title(c("3rd sample with ","sample mean ", mean),c("darkgreen","green3","green3"))
  })
 
  output$classiccltmeans<-renderPlot({        #plot of sample means based on distribution and then on whether the x-axis is the same as the 3 samples' x-axis - reference comments in the normal sample mean plot for other distributions' sample mean plot
    if (input$simd == "Normal") {
      set.seed(input$seed) #ensures that we'll see the same means from the 3 sample plots
      storage<-rep(NA,as.numeric(input$N)) #create an object to store the means of however many samples the user inputted
      for(ii in 1:as.numeric(input$N)) {
        storage[ii]<-mean(rnorm(as.numeric(input$n),as.numeric(input$norm1),as.numeric(input$norm2))) #get the means
      }
      if (input$axis == "Yes") { #if this plot's x-axis is the same as the x-axis of the sample plots then: 
          data1 <- first_data()  # Access the value of 'first' from the reactiveVal
          data2 <- second_data()  
          data3 <- third_data()  
          xmin<-min(data1,data2,data3)
          xmax<-max(data1,data2,data3)
          hist<-hist(storage,main="",
                xlab="Sample mean",col="cornflowerblue",xlim=c(floor(xmin),ceiling(xmax))) #same x-axis limits as sample plots
          multicolor.title(c("Sample means of ",input$N, " normal samples with", " normal curve"), c("cornflowerblue","cornflowerblue","cornflowerblue","black"))
          
          data1 <- first_data()  # Access means of each sample and add them on to plot in the same color
          data2 <- second_data()  
          data3 <- third_data()  
          abline(v=mean(data1),col="mediumorchid1",lwd=5)
          abline(v=mean(data2),col="salmon1",lwd=5)
          abline(v=mean(data3),col="green3",lwd=5)
          hist
          
          #make approximate normal curve
          data_mean <- mean(storage) #get mean and sd of sample means
          data_sd <- sd(storage)
          
          # Calculate x values for the normal curve
          x <- seq(floor(xmin), ceiling(xmax), length.out = 10000*as.numeric(input$n)) #the more values there are, the smoother the curve
          
          # Calculate the corresponding y values from the normal distribution
          y <- dnorm(x, mean = data_mean, sd = data_sd) 
          
          # Get the maximum count of the histogram
          max_count <- max(hist$counts)
          
          # Scale the y-values of the normal curve to match the maximum count of the histogram
          y_scaled <- y * (max_count / max(y))
          
          # Add the normal curve to the plot
          lines(x, y_scaled, col ="black", lwd = 2)
          
      } else { #if you don't want the same x-axis as the 3 sample plots then: 
        hist<-hist(storage,main="",
                   xlab="Sample mean",col="cornflowerblue") #plot normally
        multicolor.title(c("Sample means of ",input$N, " normal samples with", " normal curve"), c("cornflowerblue","cornflowerblue","cornflowerblue","black"))
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        abline(v=mean(data1),col="mediumorchid1",lwd=5)
        abline(v=mean(data2),col="salmon1",lwd=5)
        abline(v=mean(data3),col="green3",lwd=5)
        hist
        
        data_mean <- mean(storage)
        data_sd <- sd(storage)
        
        # Calculate x values for the normal curve
        x <- seq(floor(min(storage)), ceiling(max(storage)), length.out = 10000*as.numeric(input$n))
        
        # Calculate the corresponding y values from the normal distribution
        y <- dnorm(x, mean = data_mean, sd = data_sd) 
        
        # Get the maximum count of the histogram
        max_count <- max(hist$counts)
        
        # Scale the y-values of the normal curve to match the maximum count of the histogram
        y_scaled <- y * (max_count / max(y))
        
        # Add the normal curve to the plot
        lines(x, y_scaled, col ="black", lwd = 2)
      }
    }
    else if (input$simd == "Negative binomial") {
      set.seed(input$seed)
      storage<-rep(NA,as.numeric(input$N))
      for(ii in 1:as.numeric(input$N)) {
        storage[ii]<-mean(rnbinom(as.numeric(input$n),as.numeric(input$nb1),as.numeric(input$nb2)))
      }
      if (input$axis == "Yes") {
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        xmin<-min(data1,data2,data3)
        xmax<-max(data1,data2,data3)
        hist<-hist(storage,main="",
                   xlab="Sample mean",col="cornflowerblue",xlim=c(floor(xmin),ceiling(xmax)))
        multicolor.title(c("Sample means of ",input$N, " negative binomial samples with", " normal curve"), c("cornflowerblue","cornflowerblue","cornflowerblue","black"))
        
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        abline(v=mean(data1),col="mediumorchid1",lwd=5)
        abline(v=mean(data2),col="salmon1",lwd=5)
        abline(v=mean(data3),col="green3",lwd=5)
        hist
        
        data_mean <- mean(storage)
        data_sd <- sd(storage)
        
        # Calculate x values for the normal curve
        x <- seq(floor(xmin), ceiling(xmax), length.out = 10000*as.numeric(input$n))
        
        # Calculate the corresponding y values from the normal distribution
        y <- dnorm(x, mean = data_mean, sd = data_sd) 
        
        # Get the maximum count of the histogram
        max_count <- max(hist$counts)
        
        # Scale the y-values of the normal curve to match the maximum count of the histogram
        y_scaled <- y * (max_count / max(y))
        
        # Add the normal curve to the plot
        lines(x, y_scaled, col ="black", lwd = 2)
        
      } else {
        hist<-hist(storage,main="",
                   xlab="Sample mean",col="cornflowerblue")
        multicolor.title(c("Sample means of ",input$N, " negative binomial samples with", " normal curve"), c("cornflowerblue","cornflowerblue","cornflowerblue","black"))
    
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        abline(v=mean(data1),col="mediumorchid1",lwd=5)
        abline(v=mean(data2),col="salmon1",lwd=5)
        abline(v=mean(data3),col="green3",lwd=5)
        hist
        
        data_mean <- mean(storage)
        data_sd <- sd(storage)
        
        # Calculate x values for the normal curve
        x <- seq(floor(min(storage)), ceiling(max(storage)), length.out = 10000*as.numeric(input$n))
        
        # Calculate the corresponding y values from the normal distribution
        y <- dnorm(x, mean = data_mean, sd = data_sd) 
        
        # Get the maximum count of the histogram
        max_count <- max(hist$counts)
        
        # Scale the y-values of the normal curve to match the maximum count of the histogram
        y_scaled <- y * (max_count / max(y))
        
        # Add the normal curve to the plot
        lines(x, y_scaled, col ="black", lwd = 2)
      }
    }
    else if (input$simd == "Beta") {
      set.seed(input$seed)
      storage<-rep(NA,as.numeric(input$N))
      for(ii in 1:as.numeric(input$N)) {
        storage[ii]<-mean(rbeta(as.numeric(input$n),as.numeric(input$beta1),as.numeric(input$beta2)))
      }
      if (input$axis == "Yes") {
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        xmin<-min(data1,data2,data3)
        xmax<-max(data1,data2,data3)
        hist<-hist(storage,main="",
                   xlab="Sample mean",col="cornflowerblue",xlim=c(floor(xmin),ceiling(xmax)))
        multicolor.title(c("Sample means of ",input$N, " beta samples with", " normal curve"), c("cornflowerblue","cornflowerblue","cornflowerblue","black"))
        
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        abline(v=mean(data1),col="mediumorchid1",lwd=5)
        abline(v=mean(data2),col="salmon1",lwd=5)
        abline(v=mean(data3),col="green3",lwd=5)
        hist
        
        data_mean <- mean(storage)
        data_sd <- sd(storage)
        
        # Calculate x values for the normal curve
        x <- seq(floor(xmin), ceiling(xmax), length.out = 10000*as.numeric(input$n))
        
        # Calculate the corresponding y values from the normal distribution
        y <- dnorm(x, mean = data_mean, sd = data_sd) 
        
        # Get the maximum count of the histogram
        max_count <- max(hist$counts)
        
        # Scale the y-values of the normal curve to match the maximum count of the histogram
        y_scaled <- y * (max_count / max(y))
        
        # Add the normal curve to the plot
        lines(x, y_scaled, col ="black", lwd = 2)
        
      } else {
        hist<-hist(storage,main="",
                   xlab="Sample mean",col="cornflowerblue")
        multicolor.title(c("Sample means of ",input$N, " beta samples with", " normal curve"), c("cornflowerblue","cornflowerblue","cornflowerblue","black"))
        
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        abline(v=mean(data1),col="mediumorchid1",lwd=5)
        abline(v=mean(data2),col="salmon1",lwd=5)
        abline(v=mean(data3),col="green3",lwd=5)
        hist
        
        data_mean <- mean(storage)
        data_sd <- sd(storage)
        
        # Calculate x values for the normal curve
        x <- seq(floor(min(storage)), ceiling(max(storage)), length.out = 10000*as.numeric(input$n))
        
        # Calculate the corresponding y values from the normal distribution
        y <- dnorm(x, mean = data_mean, sd = data_sd) 
        
        # Get the maximum count of the histogram
        max_count <- max(hist$counts)
        
        # Scale the y-values of the normal curve to match the maximum count of the histogram
        y_scaled <- y * (max_count / max(y))
        
        # Add the normal curve to the plot
        lines(x, y_scaled, col ="black", lwd = 2)
      }
    }
    else if (input$simd == "Uniform") {
      set.seed(input$seed)
      storage<-rep(NA,as.numeric(input$N))
      for(ii in 1:as.numeric(input$N)) {
        storage[ii]<-mean(runif(as.numeric(input$n),as.numeric(input$unif1),as.numeric(input$unif2)))
      }
      if (input$axis == "Yes") {
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        xmin<-min(data1,data2,data3)
        xmax<-max(data1,data2,data3)
        hist<-hist(storage,main="",
                   xlab="Sample mean",col="cornflowerblue",xlim=c(floor(xmin),ceiling(xmax)))
        multicolor.title(c("Sample means of ",input$N, " uniform samples with", " normal curve"), c("cornflowerblue","cornflowerblue","cornflowerblue","black"))
        
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        abline(v=mean(data1),col="mediumorchid1",lwd=5)
        abline(v=mean(data2),col="salmon1",lwd=5)
        abline(v=mean(data3),col="green3",lwd=5)
        hist
        
        data_mean <- mean(storage)
        data_sd <- sd(storage)
        
        # Calculate x values for the normal curve
        x <- seq(floor(xmin), ceiling(xmax), length.out = 10000*as.numeric(input$n))
        
        # Calculate the corresponding y values from the normal distribution
        y <- dnorm(x, mean = data_mean, sd = data_sd) 
        
        # Get the maximum count of the histogram
        max_count <- max(hist$counts)
        
        # Scale the y-values of the normal curve to match the maximum count of the histogram
        y_scaled <- y * (max_count / max(y))
        
        # Add the normal curve to the plot
        lines(x, y_scaled, col ="black", lwd = 2)
        
      } else {
        hist<-hist(storage,main="",
                   xlab="Sample mean",col="cornflowerblue")
        multicolor.title(c("Sample means of ",input$N, " uniform samples with", " normal curve"), c("cornflowerblue","cornflowerblue","cornflowerblue","black"))
        
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        abline(v=mean(data1),col="mediumorchid1",lwd=5)
        abline(v=mean(data2),col="salmon1",lwd=5)
        abline(v=mean(data3),col="green3",lwd=5)
        hist
        
        data_mean <- mean(storage)
        data_sd <- sd(storage)
        
        # Calculate x values for the normal curve
        x <- seq(floor(min(storage)), ceiling(max(storage)), length.out = 10000*as.numeric(input$n))
        
        # Calculate the corresponding y values from the normal distribution
        y <- dnorm(x, mean = data_mean, sd = data_sd) 
        
        # Get the maximum count of the histogram
        max_count <- max(hist$counts)
        
        # Scale the y-values of the normal curve to match the maximum count of the histogram
        y_scaled <- y * (max_count / max(y))
        
        # Add the normal curve to the plot
        lines(x, y_scaled, col ="black", lwd = 2)
      }
      
    }
    else if (input$simd == "Poisson") {
      set.seed(input$seed)
      storage<-rep(NA,as.numeric(input$N))
      for(ii in 1:as.numeric(input$N)) {
        storage[ii]<-mean(rpois(as.numeric(input$n),as.numeric(input$pois1)))
      }
      if (input$axis == "Yes") {
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        xmin<-min(data1,data2,data3)
        xmax<-max(data1,data2,data3)
        hist<-hist(storage,main="",
                   xlab="Sample mean",col="cornflowerblue",xlim=c(floor(xmin),ceiling(xmax)))
        multicolor.title(c("Sample means of ",input$N, " Poisson samples with", " normal curve"), c("cornflowerblue","cornflowerblue","cornflowerblue","black"))
        
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        abline(v=mean(data1),col="mediumorchid1",lwd=5)
        abline(v=mean(data2),col="salmon1",lwd=5)
        abline(v=mean(data3),col="green3",lwd=5)
        hist
        
        data_mean <- mean(storage)
        data_sd <- sd(storage)
        
        # Calculate x values for the normal curve
        x <- seq(floor(xmin), ceiling(xmax), length.out = 10000*as.numeric(input$n))
        
        # Calculate the corresponding y values from the normal distribution
        y <- dnorm(x, mean = data_mean, sd = data_sd) 
        
        # Get the maximum count of the histogram
        max_count <- max(hist$counts)
        
        # Scale the y-values of the normal curve to match the maximum count of the histogram
        y_scaled <- y * (max_count / max(y))
        
        # Add the normal curve to the plot
        lines(x, y_scaled, col ="black", lwd = 2)
        
      } else {
        hist<-hist(storage,main="",
                   xlab="Sample mean",col="cornflowerblue")
        multicolor.title(c("Sample means of ",input$N, " Poisson samples with", " normal curve"), c("cornflowerblue","cornflowerblue","cornflowerblue","black"))
        
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        abline(v=mean(data1),col="mediumorchid1",lwd=5)
        abline(v=mean(data2),col="salmon1",lwd=5)
        abline(v=mean(data3),col="green3",lwd=5)
        hist
        
        data_mean <- mean(storage)
        data_sd <- sd(storage)
        
        # Calculate x values for the normal curve
        x <- seq(floor(min(storage)), ceiling(max(storage)), length.out = 10000*as.numeric(input$n))
        
        # Calculate the corresponding y values from the normal distribution
        y <- dnorm(x, mean = data_mean, sd = data_sd) 
        
        # Get the maximum count of the histogram
        max_count <- max(hist$counts)
        
        # Scale the y-values of the normal curve to match the maximum count of the histogram
        y_scaled <- y * (max_count / max(y))
        
        # Add the normal curve to the plot
        lines(x, y_scaled, col ="black", lwd = 2)
      }
    }
    else if (input$simd == "Geometric") {
      set.seed(input$seed)
      storage<-rep(NA,as.numeric(input$N))
      for(ii in 1:as.numeric(input$N)) {
        storage[ii]<-mean(rgeom(as.numeric(input$n),as.numeric(input$geom1)))
      }
      if (input$axis == "Yes") {
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        xmin<-min(data1,data2,data3)
        xmax<-max(data1,data2,data3)
        hist<-hist(storage,main="",
                   xlab="Sample mean",col="cornflowerblue",xlim=c(floor(xmin),ceiling(xmax)))
        multicolor.title(c("Sample means of ",input$N, " geometric samples with", " normal curve"), c("cornflowerblue","cornflowerblue","cornflowerblue","black"))
        
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        abline(v=mean(data1),col="mediumorchid1",lwd=5)
        abline(v=mean(data2),col="salmon1",lwd=5)
        abline(v=mean(data3),col="green3",lwd=5)
        hist
        
        data_mean <- mean(storage)
        data_sd <- sd(storage)
        
        # Calculate x values for the normal curve
        x <- seq(floor(xmin), ceiling(xmax), length.out = 10000*as.numeric(input$n))
        
        # Calculate the corresponding y values from the normal distribution
        y <- dnorm(x, mean = data_mean, sd = data_sd) 
        
        # Get the maximum count of the histogram
        max_count <- max(hist$counts)
        
        # Scale the y-values of the normal curve to match the maximum count of the histogram
        y_scaled <- y * (max_count / max(y))
        
        # Add the normal curve to the plot
        lines(x, y_scaled, col ="black", lwd = 2)
        
      } else {
        hist<-hist(storage,main="",
                   xlab="Sample mean",col="cornflowerblue")
        multicolor.title(c("Sample means of ",input$N, " geometric samples with", " normal curve"), c("cornflowerblue","cornflowerblue","cornflowerblue","black"))
        
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        abline(v=mean(data1),col="mediumorchid1",lwd=5)
        abline(v=mean(data2),col="salmon1",lwd=5)
        abline(v=mean(data3),col="green3",lwd=5)
        hist
        
        data_mean <- mean(storage)
        data_sd <- sd(storage)
        
        # Calculate x values for the normal curve
        x <- seq(floor(min(storage)), ceiling(max(storage)), length.out = 10000*as.numeric(input$n))
        
        # Calculate the corresponding y values from the normal distribution
        y <- dnorm(x, mean = data_mean, sd = data_sd) 
        
        # Get the maximum count of the histogram
        max_count <- max(hist$counts)
        
        # Scale the y-values of the normal curve to match the maximum count of the histogram
        y_scaled <- y * (max_count / max(y))
        
        # Add the normal curve to the plot
        lines(x, y_scaled, col ="black", lwd = 2)
      }
    }
    else if (input$simd == "Exponential") {
      set.seed(input$seed)
      storage<-rep(NA,as.numeric(input$N))
      for(ii in 1:as.numeric(input$N)) {
        storage[ii]<-mean(rexp(as.numeric(input$n),as.numeric(input$exp1)))
      }
      if (input$axis == "Yes") {
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        xmin<-min(data1,data2,data3)
        xmax<-max(data1,data2,data3)
        hist<-hist(storage,main="",
                   xlab="Sample mean",col="cornflowerblue",xlim=c(floor(xmin),ceiling(xmax)))
        multicolor.title(c("Sample means of ",input$N, " exponential samples with", " normal curve"), c("cornflowerblue","cornflowerblue","cornflowerblue","black"))
        
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        abline(v=mean(data1),col="mediumorchid1",lwd=5)
        abline(v=mean(data2),col="salmon1",lwd=5)
        abline(v=mean(data3),col="green3",lwd=5)
        hist
        
        data_mean <- mean(storage)
        data_sd <- sd(storage)
        
        # Calculate x values for the normal curve
        x <- seq(floor(xmin), ceiling(xmax), length.out = 10000*as.numeric(input$n))
        
        # Calculate the corresponding y values from the normal distribution
        y <- dnorm(x, mean = data_mean, sd = data_sd) 
        
        # Get the maximum count of the histogram
        max_count <- max(hist$counts)
        
        # Scale the y-values of the normal curve to match the maximum count of the histogram
        y_scaled <- y * (max_count / max(y))
        
        # Add the normal curve to the plot
        lines(x, y_scaled, col ="black", lwd = 2)
      
      } else {
        hist<-hist(storage,main="",
                   xlab="Sample mean",col="cornflowerblue")
        multicolor.title(c("Sample means of ",input$N, " exponential samples with", " normal curve"), c("cornflowerblue","cornflowerblue","cornflowerblue","black"))
        
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        abline(v=mean(data1),col="mediumorchid1",lwd=5)
        abline(v=mean(data2),col="salmon1",lwd=5)
        abline(v=mean(data3),col="green3",lwd=5)
        hist
        
        data_mean <- mean(storage)
        data_sd <- sd(storage)
        
        # Calculate x values for the normal curve
        x <- seq(floor(min(storage)), ceiling(max(storage)), length.out = 10000*as.numeric(input$n))
        
        # Calculate the corresponding y values from the normal distribution
        y <- dnorm(x, mean = data_mean, sd = data_sd) 
        
        # Get the maximum count of the histogram
        max_count <- max(hist$counts)
        
        # Scale the y-values of the normal curve to match the maximum count of the histogram
        y_scaled <- y * (max_count / max(y))
        
        # Add the normal curve to the plot
        lines(x, y_scaled, col ="black", lwd = 2)
      }
    }
    else if (input$simd == "Hypergeometric") {
      set.seed(input$seed)
      storage<-rep(NA,as.numeric(input$N))
      for(ii in 1:as.numeric(input$N)) {
        storage[ii]<-mean(rhyper(as.numeric(input$n),as.numeric(input$hyp1),as.numeric(input$hyp2),as.numeric(input$hyp3)))
      }
      if (input$axis == "Yes") {
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        xmin<-min(data1,data2,data3)
        xmax<-max(data1,data2,data3)
        hist<-hist(storage,main="",
                   xlab="Sample mean",col="cornflowerblue",xlim=c(floor(xmin),ceiling(xmax)))
        multicolor.title(c("Sample means of ",input$N, " hypergeometric samples with", " normal curve"), c("cornflowerblue","cornflowerblue","cornflowerblue","black"))
        
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        abline(v=mean(data1),col="mediumorchid1",lwd=5)
        abline(v=mean(data2),col="salmon1",lwd=5)
        abline(v=mean(data3),col="green3",lwd=5)
        hist
        
        data_mean <- mean(storage)
        data_sd <- sd(storage)
        
        # Calculate x values for the normal curve
        x <- seq(floor(xmin), ceiling(xmax), length.out = 10000*as.numeric(input$n))
        
        # Calculate the corresponding y values from the normal distribution
        y <- dnorm(x, mean = data_mean, sd = data_sd) 
        
        # Get the maximum count of the histogram
        max_count <- max(hist$counts)
        
        # Scale the y-values of the normal curve to match the maximum count of the histogram
        y_scaled <- y * (max_count / max(y))
        
        # Add the normal curve to the plot
        lines(x, y_scaled, col ="black", lwd = 2)
        
      } else {
        hist<-hist(storage,main="",
                   xlab="Sample mean",col="cornflowerblue")
        multicolor.title(c("Sample means of ",input$N, " hypergeometric samples with", " normal curve"), c("cornflowerblue","cornflowerblue","cornflowerblue","black"))
        
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        abline(v=mean(data1),col="mediumorchid1",lwd=5)
        abline(v=mean(data2),col="salmon1",lwd=5)
        abline(v=mean(data3),col="green3",lwd=5)
        hist
        
        data_mean <- mean(storage)
        data_sd <- sd(storage)
        
        # Calculate x values for the normal curve
        x <- seq(floor(min(storage)), ceiling(max(storage)), length.out = 10000*as.numeric(input$n))
        
        # Calculate the corresponding y values from the normal distribution
        y <- dnorm(x, mean = data_mean, sd = data_sd) 
        
        # Get the maximum count of the histogram
        max_count <- max(hist$counts)
        
        # Scale the y-values of the normal curve to match the maximum count of the histogram
        y_scaled <- y * (max_count / max(y))
        
        # Add the normal curve to the plot
        lines(x, y_scaled, col ="black", lwd = 2)
      }
    }
    else if (input$simd == "Binomial") {
      set.seed(input$seed)
      storage<-rep(NA,as.numeric(input$N))
      for(ii in 1:as.numeric(input$N)) {
        storage[ii]<-mean(rbinom(as.numeric(input$n),as.numeric(input$bin1),as.numeric(input$bin2)))
      }
      if (input$axis == "Yes") {
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        xmin<-min(data1,data2,data3)
        xmax<-max(data1,data2,data3)
        hist<-hist(storage,main="",
                   xlab="Sample mean",col="cornflowerblue",xlim=c(floor(xmin),ceiling(xmax)))
        multicolor.title(c("Sample means of ",input$N, " binomial samples with", " normal curve"), c("cornflowerblue","cornflowerblue","cornflowerblue","black"))
        
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        abline(v=mean(data1),col="mediumorchid1",lwd=5)
        abline(v=mean(data2),col="salmon1",lwd=5)
        abline(v=mean(data3),col="green3",lwd=5)
        hist
        
        data_mean <- mean(storage)
        data_sd <- sd(storage)
        
        # Calculate x values for the normal curve
        x <- seq(floor(xmin), ceiling(xmax), length.out = 10000*as.numeric(input$n))
        
        # Calculate the corresponding y values from the normal distribution
        y <- dnorm(x, mean = data_mean, sd = data_sd) 
        
        # Get the maximum count of the histogram
        max_count <- max(hist$counts)
        
        # Scale the y-values of the normal curve to match the maximum count of the histogram
        y_scaled <- y * (max_count / max(y))
        
        # Add the normal curve to the plot
        lines(x, y_scaled, col ="black", lwd = 2)
        
      } else {
        hist<-hist(storage,main="",
                   xlab="Sample mean",col="cornflowerblue")
        multicolor.title(c("Sample means of ",input$N, " binomial samples with", " normal curve"), c("cornflowerblue","cornflowerblue","cornflowerblue","black"))
        
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        abline(v=mean(data1),col="mediumorchid1",lwd=5)
        abline(v=mean(data2),col="salmon1",lwd=5)
        abline(v=mean(data3),col="green3",lwd=5)
        hist
        
        data_mean <- mean(storage)
        data_sd <- sd(storage)
        
        # Calculate x values for the normal curve
        x <- seq(floor(min(storage)), ceiling(max(storage)), length.out = 10000*as.numeric(input$n))
        
        # Calculate the corresponding y values from the normal distribution
        y <- dnorm(x, mean = data_mean, sd = data_sd) 
        
        # Get the maximum count of the histogram
        max_count <- max(hist$counts)
        
        # Scale the y-values of the normal curve to match the maximum count of the histogram
        y_scaled <- y * (max_count / max(y))
        
        # Add the normal curve to the plot
        lines(x, y_scaled, col ="black", lwd = 2)
      }
    }
    else {
      set.seed(input$seed)
      storage<-rep(NA,as.numeric(input$N))
      for(ii in 1:as.numeric(input$N)) {
        storage[ii]<-mean(rcauchy(as.numeric(input$n),as.numeric(input$cauchy1),as.numeric(input$cauchy2)))
      }
      if (input$axis == "Yes") {
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        xmin<-min(data1,data2,data3)
        xmax<-max(data1,data2,data3)
        hist<-hist(storage,main="",
                   xlab="Sample mean",col="cornflowerblue",xlim=c(floor(xmin),ceiling(xmax)))
        multicolor.title(c("Sample means of ",input$N, " Cauchy samples with", " normal curve"), c("cornflowerblue","cornflowerblue","cornflowerblue","black"))
        
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        abline(v=mean(data1),col="mediumorchid1",lwd=5)
        abline(v=mean(data2),col="salmon1",lwd=5)
        abline(v=mean(data3),col="green3",lwd=5)
        hist
        
        data_mean <- mean(storage)
        data_sd <- sd(storage)
        
        # Calculate x values for the normal curve
        x <- seq(floor(xmin), ceiling(xmax), length.out = 10000*as.numeric(input$n))
        
        # Calculate the corresponding y values from the normal distribution
        y <- dnorm(x, mean = data_mean, sd = data_sd) 
        
        # Get the maximum count of the histogram
        max_count <- max(hist$counts)
        
        # Scale the y-values of the normal curve to match the maximum count of the histogram
        y_scaled <- y * (max_count / max(y))
        
        # Add the normal curve to the plot
        lines(x, y_scaled, col ="black", lwd = 2)
        
      } else {
        hist<-hist(storage,main="",
                   xlab="Sample mean",col="cornflowerblue")
        multicolor.title(c("Sample means of ",input$N, " Cauchy samples with", " normal curve"), c("cornflowerblue","cornflowerblue","cornflowerblue","black"))
        
        data1 <- first_data()  # Access the value of 'first' from the reactiveVal
        data2 <- second_data()  
        data3 <- third_data()  
        abline(v=mean(data1),col="mediumorchid1",lwd=5)
        abline(v=mean(data2),col="salmon1",lwd=5)
        abline(v=mean(data3),col="green3",lwd=5)
        hist
        
        data_mean <- mean(storage)
        data_sd <- sd(storage)
        
        # Calculate x values for the normal curve
        x <- seq(hist$breaks[1], hist$breaks[length(hist$breaks)], length.out = 10000*as.numeric(input$n)) #hist$breaks covers the entire x-axis unlike floor(min(storage)) and ceiling(max(storage))
        
        # Calculate the corresponding y values from the normal distribution
        y <- dnorm(x, mean = data_mean, sd = data_sd) 
        
        # Get the maximum count of the histogram
        max_count <- max(hist$counts)
        
        # Scale the y-values of the normal curve to match the maximum count of the histogram
        y_scaled <- y * (max_count / max(y))
        
        # Add the normal curve to the plot
        lines(x, y_scaled, col ="black", lwd = 2)
      }
    }
  })
  output$select1 <- renderUI({                #choose 1st parameter based on chosen distribution for K-S simulation
      if (input$cltd == "Normal") {
      textInput("anorm1", "Enter mean","0")
      }
      else if (input$cltd == "Negative binomial") {
      textInput("anb1", "Enter number of successes (>0)","3")
      }
      else if (input$cltd == "Beta") {
      textInput("abeta1", "Enter α (>0)",".5")
      }
      else if (input$cltd == "Uniform") {
      textInput("aunif1", "Enter a","0")
      }
      else if (input$cltd == "Poisson") {
      textInput("apois1", "Enter λ (>0)","1")
      }
      else if (input$cltd == "Geometric") {
      textInput("ageom1", "Enter probability of success [0,1]",".5")
      }
      else if (input$cltd == "Exponential") {
      textInput("aexp1", "Enter λ (>0)","1")
      }
      else if (input$cltd == "Hypergeometric") {
      textInput("ahyp1", "Enter number of white balls (w) in urn","1")
      }
      else if (input$cltd == "Binomial") {
      textInput("abin1", "Enter number of trials (integer ≥ 0)","10")
      }
      else {
        textInput("acauchy1", "Enter location","0")
      }
  })
  output$select2 <- renderUI({             #choose 2nd parameter based on chosen distribution for K-S simulation
    if (input$cltd == "Normal") {
      textInput("anorm2", "Enter standard deviation","1")
    }
    else if (input$cltd == "Negative binomial") {
      textInput("anb2", "Enter probability of success [0,1]",".5")
    }
    else if (input$cltd == "Beta") {
      textInput("abeta2", "Enter β (>0)",".5")
    }
    else if (input$cltd == "Uniform") {
      textInput("aunif2", "Enter b (>a)","1")
    }
    else if (input$cltd == "Hypergeometric") {
      textInput("ahyp2", "Enter number of black balls (b) in urn","1")
    }
    else if (input$cltd == "Binomial") {
      textInput("abin2", "Enter probability of success [0,1]",".5")
    }
    else if (input$cltd == "Cauchy"){
      textInput("acauchy2", "Enter scale","1")
    }
    else {
    }
  })
  output$select3 <- renderUI({                   #choose 3rd parameter based on chosen distribution for K-S simulation
    if (input$cltd == "Hypergeometric") {
      textInput("ahyp3", "Enter number of balls drawn from urn (w+b ≥ integer ≥ 0)","1")
    }
    else {
    }
  })

  table_data <- reactiveValues(data = data.frame(matrix(nrow = 6, ncol = 0))) #create reactive table- starts with 0 columns
  
  col_counter <- reactiveValues(index = 1) # create reactive column counter with index 1
  max_columns <- 7  # Set your desired maximum number of columns
  
  
  # Add event handler for the "Update table and graph" button
  observeEvent(input$update, {
    if (col_counter$index <= max_columns) {   #if there are <= 7 columns then continue adding columns
    
    # Create a new column name based on the inputs that describe the distribution and its parameters
    col_name <- if (input$cltd == "Normal") {
      paste("N(", input$anorm1,",",input$anorm2,")")
    }
    else if (input$cltd == "Negative binomial") {
      paste("NBinom(", input$anb1,",",input$anb2,")")
    }
    else if (input$cltd == "Beta") {
      paste("Beta(", input$abeta1,",",input$abeta2,")")
    }
    else if (input$cltd == "Uniform") {
      paste("Unif(", input$aunif1,",",input$aunif2,")")
    }
    else if (input$cltd == "Poisson") {
      paste("Pois(", input$apois1,")")
    }
    else if (input$cltd == "Geometric") {
      paste("Geom(", input$ageom1,")")
    }
    else if (input$cltd == "Exponential") {
      paste("Expo(", input$aexp1,")")
    }
    else if (input$cltd == "Hypergeometric") {
      paste("Hyper(", input$ahyp1,",",input$ahyp2,",",input$ahyp3,")")
    }
    else if (input$cltd == "Binomial") {
      paste("Binom(", input$abin1,",",input$abin2,")")
    }
    else {
      paste("Cauchy(", input$acauchy1,",",input$acauchy2,")")
    }
    
    # Update the table data by adding a new column - reference comments for the Normal distribution column for explanation on other distributions' columns
    table_data$data[, col_name] <- if (input$cltd == "Normal") {
      #theoretical K-S
      n<-500
      m<-as.numeric(input$anorm1)
      s<-as.numeric(input$anorm2)
      norm<-rnorm(n,m,s)
      
      seq <- seq(min(norm), max(norm), length.out = n)
      cdf<-pnorm(seq,m,s)
      normcdf<-pnorm(seq,m,s)
      
      diff<-abs(cdf-normcdf)
      th.ks<-max(diff)
      th.ks<-round(th.ks,3)
      
      #sample K-S statistics - create function that gives sample K-S statistic based on sample size
      normal.ks = function(nn)
      {
        storage<-rep(NA,10000) #store 10,000 means
        for (ii in 1:10000) {
          storage[ii]<-mean(rnorm(nn,m,s))
        }
        ks = ks.test(storage, pnorm, m, sqrt(s^2/nn))$statistic #sample variance needs to be divided by sample size
        ks<-round(ks,3)
        return(ks) #KS statistic
      }
      ks0<-normal.ks(1)
      ks1<-normal.ks(10)
      ks2<-normal.ks(100)
      ks3<-normal.ks(1000)
      ks4<-normal.ks(10000)
      
      #column returns theoretical KS statistic and KS statistics for n=1,10,100,1000,10000
      c(th.ks, ks0, ks1, ks2, ks3, ks4)
    }
    else if (input$cltd == "Negative binomial") {
      #theoretical K-S
      n<-500
      r<-as.numeric(input$anb1)
      p<-as.numeric(input$anb2)
      nbinom<-rnbinom(n,r,p)
      
      seq <- seq(min(nbinom), max(nbinom), length.out = n)
      cdf<-pnbinom(seq,r,p)
      normcdf<-pnorm(seq,(r*(1-p))/p,sqrt((r*(1-p))/(p^2)))
      
      diff<-abs(cdf-normcdf)
      th.ks<-max(diff) 
      th.ks<-round(th.ks,3)
      
      #sample K-S
      fun.ks = function(nn)
      {
        storage<-rep(NA,10000)
        for (ii in 1:10000) {
          storage[ii]<-mean(rnbinom(nn,r,p))
        }
        ks = ks.test(storage, pnorm, (r*(1-p))/p, sqrt(((r*(1-p))/(p^2))/nn))$statistic
        ks<-round(ks,3)
        return(ks)
      }
      ks0<-fun.ks(1)
      ks1<-fun.ks(10)
      ks2<-fun.ks(100)
      ks3<-fun.ks(1000)
      ks4<-fun.ks(10000)
      
      c(th.ks, ks0,ks1,ks2,ks3,ks4)
    }
    else if (input$cltd == "Beta") {
      n<-500
      a<-as.numeric(input$abeta1)
      b<-as.numeric(input$abeta2)
      beta<-rbeta(n,a,b)
      
      seq <- seq(min(beta), max(beta), length.out = n)
      cdf<-pbeta(seq,a,b)
      normcdf<-pnorm(seq,a/(a+b),sqrt((a*b)/((a+b)^2*(a+b+1))))
      
      diff<-abs(cdf-normcdf)
      th.ks<-max(diff) #theoretical K-S stat
      th.ks<-round(th.ks,3)
      
      beta.ks = function(nn)
      {
        storage<-rep(NA,10000) 
        for (ii in 1:10000) {
          storage[ii]<-mean(rbeta(nn, a, b)) 
        }
        ks=ks.test(storage, pnorm, a/(a+b),sqrt(((a*b)/((a+b)^2*(a+b+1)))/nn))$statistic
        return(ks)
      }
      ks0<-beta.ks(1)
      ks1<-beta.ks(10)
      ks2<-beta.ks(100)
      ks3<-beta.ks(1000)
      ks4<-beta.ks(10000)
      
      c(th.ks, ks0,ks1,ks2,ks3,ks4)
      
    }
    else if (input$cltd == "Uniform") {
      n<-500
      a<-as.numeric(input$aunif1)
      b<-as.numeric(input$aunif2)
      uni<-runif(n,a,b)
      
      seq <- seq(min(uni), max(uni), length.out = n)
      cdf<-punif(seq,a,b)
      normcdf<-pnorm(seq,.5*(a+b),sqrt((1/12)*(b-a)^2))
      
      diff<-abs(cdf-normcdf)
      th.ks<-max(diff) #theoretical K-S stat
      th.ks<-round(th.ks,3)
      
      fun.ks = function(nn) {
        storage<-rep(NA,10000) 
        for (ii in 1:10000) {
          storage[ii]<-mean(runif(nn,a,b))
        }
        ks=ks.test(storage, pnorm, .5*(a+b),sqrt(((1/12)*(b-a)^2)/nn))$statistic
        return(ks)
      }
      
      ks0<-fun.ks(1)
      ks1<-fun.ks(10)
      ks2<-fun.ks(100)
      ks3<-fun.ks(1000)
      ks4<-fun.ks(10000)
      
      c(th.ks, ks0,ks1,ks2,ks3,ks4)
      
    }
    else if (input$cltd == "Poisson") {
      n<-500
      l<-as.numeric(input$apois1)
      pois<-rpois(n,l)
      
      seq <- seq(min(pois), max(pois), length.out = n)
      cdf<-ppois(seq,l)
      normcdf<-pnorm(seq,l,sqrt(l))
      
      diff<-abs(cdf-normcdf)
      th.ks<-max(diff) #theoretical K-S stat
      th.ks<-round(th.ks,3)
      
      pois.ks = function(nn) {
        storage<-rep(NA,10000) 
        for (ii in 1:10000) {
          storage[ii]<-mean(rpois(nn, l)) 
        }
        ks=ks.test(storage, pnorm, l, sqrt(l/nn))$statistic
        return(ks)
      }
      
      ks0<-pois.ks(1)
      ks1<-pois.ks(10)
      ks2<-pois.ks(100)
      ks3<-pois.ks(1000)
      ks4<-pois.ks(10000)
      
      c(th.ks, ks0,ks1,ks2,ks3,ks4)
    }
    else if (input$cltd == "Geometric") {
      n<-500
      p<-as.numeric(input$ageom1)
      geom<-rgeom(n,p)
      
      seq <- seq(min(geom), max(geom), length.out = n)
      cdf<-pgeom(seq,p)
      normcdf<-pnorm(seq,(1-p)/p, sqrt((1-p)/p^2))
      
      diff<-abs(cdf-normcdf)
      th.ks<-max(diff) #theoretical K-S stat
      th.ks<-round(th.ks,3)
      
      geom.ks = function(nn)
      {
        storage<-rep(NA,10000) 
        for (ii in 1:10000) {
          storage[ii]<-mean(rgeom(nn,p)) 
        }
        ks=ks.test(storage, pnorm, (1-p)/p, sqrt(((1-p)/p^2)/nn))$statistic
        return(ks)
      }
      
      ks0<-geom.ks(1)
      ks1<-geom.ks(10)
      ks2<-geom.ks(100)
      ks3<-geom.ks(1000)
      ks4<-geom.ks(10000)
      
      c(th.ks, ks0,ks1,ks2,ks3,ks4)
    }
    else if (input$cltd == "Exponential") {
      n<-500
      l<-as.numeric(input$aexp1)
      exp<-rexp(n,l)
      
      seq <- seq(min(exp), max(exp), length.out = n)
      cdf<-pexp(seq,l)
      normcdf<-pnorm(seq,1/l, 1/l)
      
      diff<-abs(cdf-normcdf)
      th.ks<-max(diff) #theoretical K-S stat
      th.ks<-round(th.ks,3)
      
      exp.ks = function(nn)
      {
        storage<-rep(NA,10000) 
        for (ii in 1:10000) {
          storage[ii]<-mean(rexp(nn,l)) 
        }
        ks=ks.test(storage, pnorm, 1/l, (1/l)/sqrt(nn))$statistic
        return(ks)
      }
      
      ks0<-exp.ks(1)
      ks1<-exp.ks(10)
      ks2<-exp.ks(100)
      ks3<-exp.ks(1000)
      ks4<-exp.ks(10000)
      
      c(th.ks, ks0,ks1,ks2,ks3,ks4)
    }
    else if (input$cltd == "Hypergeometric") {
      n<-500
      w<-as.numeric(input$ahyp1)
      b<-as.numeric(input$ahyp2)
      t<-as.numeric(input$ahyp3)
      hyp<-rhyper(n,w,b,t)
      
      seq <- seq(min(hyp), max(hyp), length.out = n)
      cdf<-phyper(seq,w,b,t)
      normcdf<-pnorm(seq,t*(w/(w+b)), sqrt(t*(w/(w+b))*(1-(w/(w+b)))*((w+b-t)/(w+b-1))))
      
      diff<-abs(cdf-normcdf)
      th.ks<-max(diff) #theoretical K-S stat
      th.ks<-round(th.ks,3)
      
      hg.ks = function(nn)
      {
        storage<-rep(NA,10000) 
        for (ii in 1:10000) {
          storage[ii]<-mean(rhyper(nn,w,b,t)) 
        }
        ks=ks.test(storage, pnorm, t*(w/(w+b)), sqrt((t*(w/(w+b))*(1-(w/(w+b)))*((w+b-t)/(w+b-1)))/nn))$statistic
        return(ks)
      }
      ks0<-hg.ks(1)
      ks1<-hg.ks(10)
      ks2<-hg.ks(100)
      ks3<-hg.ks(1000)
      ks4<-hg.ks(10000)
      
      c(th.ks, ks0,ks1,ks2,ks3,ks4)
    }
    else if (input$cltd == "Binomial") {
      n<-500
      t<-as.numeric(input$abin1)
      p<-as.numeric(input$abin2)
      bin<-rbinom(n,t,p)
      
      seq <- seq(min(bin), max(bin), length.out = n)
      cdf<-pbinom(seq,t,p)
      normcdf<-pnorm(seq,t*p, sqrt(t*p*(1-p)))
      
      diff<-abs(cdf-normcdf)
      th.ks<-max(diff) #theoretical K-S stat
      th.ks<-round(th.ks,3)
      
      bin.ks = function(nn)
      {
        storage<-rep(NA,10000) 
        for (ii in 1:10000) {
          storage[ii]<-mean(rbinom(nn,t,p))
        }
        ks = ks.test(storage, pnorm, t*p, sqrt((t*p*(1-p))/nn))$statistic
        return(ks)
      }
      
      ks0<-bin.ks(1)
      ks1<-bin.ks(10)
      ks2<-bin.ks(100)
      ks3<-bin.ks(1000)
      ks4<-bin.ks(10000)
      
      c(th.ks, ks0,ks1,ks2,ks3,ks4)
    }
    else {
      n<-500
      l<-as.numeric(input$acauchy1)
      s<-as.numeric(input$acauchy2)
      cauchy<-rcauchy(n,l,s)
      
      #See comment on Cauchy line 717
      # sims<-1000
      # 
      # storagem<-rep(NA,sims) 
      # storagev<-rep(NA,sims)
      # for (ii in 1:sims) {
      #   storagem[ii]<-mean(rcauchy(n))
      #   storagev[ii]<-var(rcauchy(n))
      #   mean<-mean(storagem)
      #   var<-mean(storagev)
      # }
      
      seq <- seq(min(cauchy), max(cauchy), length.out = n)
      cdf<-pcauchy(seq,l,s)
      normcdf<-pnorm(seq,l, s)
      
      diff<-abs(cdf-normcdf)
      th.ks<-max(diff) #theoretical K-S stat
      th.ks<-round(th.ks,3)
      
      c.ks = function(nn)
      {
        storage<-rep(NA,10000) 
        for (ii in 1:10000) {
          storage[ii]<-mean(rcauchy(nn,l,s))
        }
        ks=ks.test(storage, pnorm, l, sqrt(s^2/nn))$statistic
        return(ks)
      }
      ks0<-c.ks(1)
      ks1<-c.ks(10)
      ks2<-c.ks(100)
      ks3<-c.ks(1000)
      ks4<-c.ks(10000)
      
      c(th.ks, ks0,ks1,ks2,ks3,ks4)
    }
  
    col_counter$index <- col_counter$index + 1 #increase column counter by 1 to represent next column
    
    } else { #but if the column added is the 8th column: 
      table_data$data <- data.frame(matrix(nrow = 6, ncol = 0)) #reset table so the new column is the first column
      col_counter$index <- 1 #reset column counter and repeat previous code
      # Create a new column name based on the inputs
      col_name <- if (input$cltd == "Normal") {
        paste("N(", input$anorm1,",",input$anorm2,")")
      }
      else if (input$cltd == "Negative binomial") {
        paste("NBinom(", input$anb1,",",input$anb2,")")
      }
      else if (input$cltd == "Beta") {
        paste("Beta(", input$abeta1,",",input$abeta2,")")
      }
      else if (input$cltd == "Uniform") {
        paste("Unif(", input$aunif1,",",input$aunif2,")")
      }
      else if (input$cltd == "Poisson") {
        paste("Pois(", input$apois1,")")
      }
      else if (input$cltd == "Geometric") {
        paste("Geom(", input$ageom1,")")
      }
      else if (input$cltd == "Exponential") {
        paste("Expo(", input$aexp1,")")
      }
      else if (input$cltd == "Hypergeometric") {
        paste("Hyper(", input$ahyp1,",",input$ahyp2,",",input$ahyp3,")")
      }
      else if (input$cltd == "Binomial") {
        paste("Binom(", input$abin1,",",input$abin2,")")
      }
      else {
        paste("Cauchy(", input$acauchy1,",",input$acauchy2,")")
      }
      
      # Update the table data by adding a new column
      table_data$data[, col_name] <- if (input$cltd == "Normal") {
        #theoretical K-S
        n<-500
        m<-as.numeric(input$anorm1)
        s<-as.numeric(input$anorm2)
        norm<-rnorm(n,m,s)
        
        seq <- seq(min(norm), max(norm), length.out = n)
        cdf<-pnorm(seq,m,s)
        normcdf<-pnorm(seq,m,s)
        
        diff<-abs(cdf-normcdf)
        th.ks<-max(diff)
        th.ks<-round(th.ks,3)
        
        #sample K-S
        normal.ks = function(nn)
        {
          storage<-rep(NA,10000)
          for (ii in 1:10000) {
            storage[ii]<-mean(rnorm(nn,m,s))
          }
          ks = ks.test(storage, pnorm, m, sqrt(s^2/nn))$statistic
          ks<-round(ks,3)
          return(ks)
        }
        ks0<-normal.ks(1)
        ks1<-normal.ks(10)
        ks2<-normal.ks(100)
        ks3<-normal.ks(1000)
        ks4<-normal.ks(10000)
        
        c(th.ks, ks0, ks1, ks2, ks3, ks4)
      }
      else if (input$cltd == "Negative binomial") {
        #theoretical K-S
        n<-500
        r<-as.numeric(input$anb1)
        p<-as.numeric(input$anb2)
        nbinom<-rnbinom(n,r,p)
        
        seq <- seq(min(nbinom), max(nbinom), length.out = n)
        cdf<-pnbinom(seq,r,p)
        normcdf<-pnorm(seq,(r*(1-p))/p,sqrt((r*(1-p))/(p^2)))
        
        diff<-abs(cdf-normcdf)
        th.ks<-max(diff) 
        th.ks<-round(th.ks,3)
        
        #sample K-S
        fun.ks = function(nn)
        {
          storage<-rep(NA,10000)
          for (ii in 1:10000) {
            storage[ii]<-mean(rnbinom(nn,r,p))
          }
          ks = ks.test(storage, pnorm, (r*(1-p))/p, sqrt(((r*(1-p))/(p^2))/nn))$statistic
          ks<-round(ks,3)
          return(ks)
        }
        ks0<-fun.ks(1)
        ks1<-fun.ks(10)
        ks2<-fun.ks(100)
        ks3<-fun.ks(1000)
        ks4<-fun.ks(10000)
        
        c(th.ks, ks0,ks1,ks2,ks3,ks4)
      }
      else if (input$cltd == "Beta") {
        n<-500
        a<-as.numeric(input$abeta1)
        b<-as.numeric(input$abeta2)
        beta<-rbeta(n,a,b)
        
        seq <- seq(min(beta), max(beta), length.out = n)
        cdf<-pbeta(seq,a,b)
        normcdf<-pnorm(seq,a/(a+b),sqrt((a*b)/((a+b)^2*(a+b+1))))
        
        diff<-abs(cdf-normcdf)
        th.ks<-max(diff) #theoretical K-S stat
        th.ks<-round(th.ks,3)
        
        beta.ks = function(nn)
        {
          storage<-rep(NA,10000) 
          for (ii in 1:10000) {
            storage[ii]<-mean(rbeta(nn, a, b)) 
          }
          ks=ks.test(storage, pnorm, a/(a+b),sqrt(((a*b)/((a+b)^2*(a+b+1)))/nn))$statistic
          return(ks)
        }
        ks0<-beta.ks(1)
        ks1<-beta.ks(10)
        ks2<-beta.ks(100)
        ks3<-beta.ks(1000)
        ks4<-beta.ks(10000)
        
        c(th.ks, ks0,ks1,ks2,ks3,ks4)
        
      }
      else if (input$cltd == "Uniform") {
        n<-500
        a<-as.numeric(input$aunif1)
        b<-as.numeric(input$aunif2)
        uni<-runif(n,a,b)
        
        seq <- seq(min(uni), max(uni), length.out = n)
        cdf<-punif(seq,a,b)
        normcdf<-pnorm(seq,.5*(a+b),sqrt((1/12)*(b-a)^2))
        
        diff<-abs(cdf-normcdf)
        th.ks<-max(diff) #theoretical K-S stat
        th.ks<-round(th.ks,3)
        
        fun.ks = function(nn) {
          storage<-rep(NA,10000) 
          for (ii in 1:10000) {
            storage[ii]<-mean(runif(nn,a,b))
          }
          ks=ks.test(storage, pnorm, .5*(a+b),sqrt(((1/12)*(b-a)^2)/nn))$statistic
          return(ks)
        }
        
        ks0<-fun.ks(1)
        ks1<-fun.ks(10)
        ks2<-fun.ks(100)
        ks3<-fun.ks(1000)
        ks4<-fun.ks(10000)
        
        c(th.ks, ks0,ks1,ks2,ks3,ks4)
        
      }
      else if (input$cltd == "Poisson") {
        n<-500
        l<-as.numeric(input$apois1)
        pois<-rpois(n,l)
        
        seq <- seq(min(pois), max(pois), length.out = n)
        cdf<-ppois(seq,l)
        normcdf<-pnorm(seq,l,sqrt(l))
        
        diff<-abs(cdf-normcdf)
        th.ks<-max(diff) #theoretical K-S stat
        th.ks<-round(th.ks,3)
        
        pois.ks = function(nn) {
          storage<-rep(NA,10000) 
          for (ii in 1:10000) {
            storage[ii]<-mean(rpois(nn, l)) 
          }
          ks=ks.test(storage, pnorm, l, sqrt(l/nn))$statistic
          return(ks)
        }
        
        ks0<-pois.ks(1)
        ks1<-pois.ks(10)
        ks2<-pois.ks(100)
        ks3<-pois.ks(1000)
        ks4<-pois.ks(10000)
        
        c(th.ks, ks0,ks1,ks2,ks3,ks4)
      }
      else if (input$cltd == "Geometric") {
        n<-500
        p<-as.numeric(input$ageom1)
        geom<-rgeom(n,p)
        
        seq <- seq(min(geom), max(geom), length.out = n)
        cdf<-pgeom(seq,p)
        normcdf<-pnorm(seq,(1-p)/p, sqrt((1-p)/p^2))
        
        diff<-abs(cdf-normcdf)
        th.ks<-max(diff) #theoretical K-S stat
        th.ks<-round(th.ks,3)
        
        geom.ks = function(nn)
        {
          storage<-rep(NA,10000) 
          for (ii in 1:10000) {
            storage[ii]<-mean(rgeom(nn,p)) 
          }
          ks=ks.test(storage, pnorm, (1-p)/p, sqrt(((1-p)/p^2)/nn))$statistic
          return(ks)
        }
        
        ks0<-geom.ks(1)
        ks1<-geom.ks(10)
        ks2<-geom.ks(100)
        ks3<-geom.ks(1000)
        ks4<-geom.ks(10000)
        
        c(th.ks, ks0,ks1,ks2,ks3,ks4)
      }
      else if (input$cltd == "Exponential") {
        n<-500
        l<-as.numeric(input$aexp1)
        exp<-rexp(n,l)
        
        seq <- seq(min(exp), max(exp), length.out = n)
        cdf<-pexp(seq,l)
        normcdf<-pnorm(seq,1/l, 1/l)
        
        diff<-abs(cdf-normcdf)
        th.ks<-max(diff) #theoretical K-S stat
        th.ks<-round(th.ks,3)
        
        exp.ks = function(nn)
        {
          storage<-rep(NA,10000) 
          for (ii in 1:10000) {
            storage[ii]<-mean(rexp(nn,l)) 
          }
          ks=ks.test(storage, pnorm, 1/l, (1/l)/sqrt(nn))$statistic
          return(ks)
        }
        
        ks0<-exp.ks(1)
        ks1<-exp.ks(10)
        ks2<-exp.ks(100)
        ks3<-exp.ks(1000)
        ks4<-exp.ks(10000)
        
        c(th.ks, ks0,ks1,ks2,ks3,ks4)
      }
      else if (input$cltd == "Hypergeometric") {
        n<-500
        w<-as.numeric(input$ahyp1)
        b<-as.numeric(input$ahyp2)
        t<-as.numeric(input$ahyp3)
        hyp<-rhyper(n,w,b,t)
        
        seq <- seq(min(hyp), max(hyp), length.out = n)
        cdf<-phyper(seq,w,b,t)
        normcdf<-pnorm(seq,t*(w/(w+b)), sqrt(t*(w/(w+b))*(1-(w/(w+b)))*((w+b-t)/(w+b-1))))
        
        diff<-abs(cdf-normcdf)
        th.ks<-max(diff) #theoretical K-S stat
        th.ks<-round(th.ks,3)
        
        hg.ks = function(nn)
        {
          storage<-rep(NA,10000) 
          for (ii in 1:10000) {
            storage[ii]<-mean(rhyper(nn,w,b,t)) 
          }
          ks=ks.test(storage, pnorm, t*(w/(w+b)), sqrt((t*(w/(w+b))*(1-(w/(w+b)))*((w+b-t)/(w+b-1)))/nn))$statistic
          return(ks)
        }
        ks0<-hg.ks(1)
        ks1<-hg.ks(10)
        ks2<-hg.ks(100)
        ks3<-hg.ks(1000)
        ks4<-hg.ks(10000)
        
        c(th.ks, ks0,ks1,ks2,ks3,ks4)
      }
      else if (input$cltd == "Binomial") {
        n<-500
        t<-as.numeric(input$abin1)
        p<-as.numeric(input$abin2)
        bin<-rbinom(n,t,p)
        
        seq <- seq(min(bin), max(bin), length.out = n)
        cdf<-pbinom(seq,t,p)
        normcdf<-pnorm(seq,t*p, sqrt(t*p*(1-p)))
        
        diff<-abs(cdf-normcdf)
        th.ks<-max(diff) #theoretical K-S stat
        th.ks<-round(th.ks,3)
        
        bin.ks = function(nn)
        {
          storage<-rep(NA,10000) 
          for (ii in 1:10000) {
            storage[ii]<-mean(rbinom(nn,t,p))
          }
          ks = ks.test(storage, pnorm, t*p, sqrt((t*p*(1-p))/nn))$statistic
          return(ks)
        }
        
        ks0<-bin.ks(1)
        ks1<-bin.ks(10)
        ks2<-bin.ks(100)
        ks3<-bin.ks(1000)
        ks4<-bin.ks(10000)
        
        c(th.ks, ks0,ks1,ks2,ks3,ks4)
      }
      else {
        n<-500
        l<-as.numeric(input$acauchy1)
        s<-as.numeric(input$acauchy2)
        cauchy<-rcauchy(n,l,s)
        
        #See comment on Cauchy on line 717
        # sims<-1000
        # 
        # storagem<-rep(NA,sims) 
        # storagev<-rep(NA,sims)
        # for (ii in 1:sims) {
        #   storagem[ii]<-mean(rcauchy(n))
        #   storagev[ii]<-var(rcauchy(n))
        #   mean<-mean(storagem)
        #   var<-mean(storagev)
        # }
        # 
        seq <- seq(min(cauchy), max(cauchy), length.out = n)
        cdf<-pcauchy(seq,l,s)
        normcdf<-pnorm(seq,l, s)
        
        diff<-abs(cdf-normcdf)
        th.ks<-max(diff) #theoretical K-S stat
        th.ks<-round(th.ks,3)
        
        c.ks = function(nn)
        {
          storage<-rep(NA,10000) 
          for (ii in 1:10000) {
            storage[ii]<-mean(rcauchy(nn,l,s))
          }
          ks=ks.test(storage, pnorm, l, sqrt(s^2/nn))$statistic
          return(ks)
        }
        ks0<-c.ks(1)
        ks1<-c.ks(10)
        ks2<-c.ks(100)
        ks3<-c.ks(1000)
        ks4<-c.ks(10000)
        
        c(th.ks, ks0,ks1,ks2,ks3,ks4)
      }
      col_counter$index <- col_counter$index + 1 #increment column counter to represent next column
      
    }
    
    #create table of K-S values
    output$classict <- renderTable({
      table_data$data <- round(table_data$data, 3) 
      formatted_data <- format(table_data$data, nsmall = 3) #show 3 digits
      rownames(formatted_data) <- c("Theoretical K-S", "n=1", "n=10", "n=100", "n=1,000", "n=10,000") #give row names
      formatted_data      
    },rownames=T)
    
    #create plot of simulated K-S values
    output$classicplot <- renderPlot({
      req(table_data$data)  # Ensure data is available
      
      # X values (natural log-scaled)
      x_values <- c(log(1), log(10), log(100), log(1000), log(10000))
      
      y_values <- table_data$data[2:6, 1]  # 2nd to 6th rows (exclude 1st row since that is the theoretical value)
      
      colors<-c("firebrick1","orange1","gold","mediumseagreen","skyblue1","slateblue4","blueviolet" ) #color vector for each line plotted

      if (ncol(table_data$data)>1) { #if there are >1 columns in the table:
        
        #then plot the log sample size and simulated K-S values in the first column corresponding to the first color in vector colors
        plot(x_values, y_values, type = "o", col = colors[1],lwd=5, 
        xlab = "Log Sample Size", ylab = "K-S Statistic",
        main = "Simulated K-S Statistics for Multiple Sample Sizes",ylim=c(0,as.numeric(input$y)))
        
        i<-2:ncol(table_data$data) #plot lines for the ith column in table corresponding to the ith color in vector colors 
        repeat {
          lines(x_values,table_data$data[2:6,i[1]],type="o",col=colors[i[1]],lwd=5)
          i<-i[-c(1)]
          if (is.na(i[1]))
            break
        }
        #add legend with the column names with the corresponding color
        legend("topright",legend=colnames(table_data$data)[1:ncol(table_data$data)], col=colors[1:ncol(table_data$data)],pch=19)
        
      }
        else { #if there's only 1 column in the table then just plot that
          plot(x_values, y_values, type = "o", col = colors[1], lwd=5,
               xlab = "Log Sample Size", ylab = "K-S Statistic",
               main = "Simulated K-S Statistics for Multiple Sample Sizes",
               ylim=c(0,as.numeric(input$y)))
          
          legend("topright",legend=colnames(table_data$data)[1], col=colors[1],pch=19)
        }
      
      })
 
 })
 
}  
shinyApp(ui=ui, server=server)
