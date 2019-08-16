##############################################
###                                        ###
###   Adaptive Design of Experiments for   ###
###       Maximizing Information Gain      ###
###        code by: Brennan Klein          ###
###                                        ###
##############################################
# Plotting Processes 

plot_landscape <- function(history, n_samples, n_init, sobol, n_dim, expParam1_range, expParam2_range, k=1, model_num=1){
	# Does the landscape plotting, without the max/min/etc. points added on.
	datagrid <- expand.grid(A= seq(expParam1_range[1], expParam1_range[2], length.out=n_dim), 
						    PI=seq(expParam2_range[1], expParam2_range[2], length.out=n_dim))
	
	coords <- jitter_points(history[c("PI","A")], datagrid, history)
	history$A  <- coords$A
	history$PI <- coords$PI
	
	GP <- run_gp(history, datagrid)
	df <- data.frame(GP$complete_data)
	colnames(df) <- c("A", "PI", "Information", "MSE")
	df <- df[df$PI>=0.2 & df$PI<=0.8,]
	max_point <- df[which.max(df$Information),]
	coords <- as.numeric(max_point[1:2])
	max_point <- data.frame(A=coords[1], PI=coords[2])

	title <- paste("Model: ", model_num, " ... Number of Searches: ", dim(history)[1], sep="")	
	
	out <- ggplot(df, aes(x=PI,y=A, guides=T)) + 
		geom_point(aes(color=Information), alpha=0.975, size=200/n_dim, shape=15) + xlim(c(0.185,0.815)) +
		scale_color_gradientn(colors=jet.colors) + 
		geom_point(data=history, aes(x=PI, y=A, fill = "Searches"), size=1.4, shape=21, color="white", alpha=0.8) +
		geom_point(data=max_point, aes(x=PI, y=A, fill = "Maximum"), size=6, shape=24, color="white", alpha=0.8) +
		ggtitle(title) + theme_minimal() + xlab(TeX("$\\pi$")) + ylab(TeX("$A$")) + 
		scale_fill_manual(name='', values=c("Searches"='black', "Maximum"='deeppink'), guide='legend') +
		guides(fill = guide_legend(override.aes = list(shape=c(24, 21), color=c("deeppink", "black"), size=c(3,2)))) + 
		theme(text=element_text(size=12, family="Times"), 
			plot.title = element_text(family="Times", hjust=0.5, size=15), 
	        axis.title.x = element_text(size=16),
	        axis.text.x = element_text(size=14),
	        axis.text.y = element_text(size=14),
	        axis.title.y = element_text(size=16),
	        axis.ticks.x=element_blank(), 
	        axis.ticks.y=element_blank(), 
	        panel.border = element_blank(), panel.grid.major = element_blank(),
	        panel.grid.minor = element_blank(), axis.line = element_blank())
	
	return(list(out))
}

make_filename <- function(asymmetric, sampling, n_init, method, n_samples, search, seed, sobol, models_used, top_model, dtype="CSV"){
	asym <- "asym"; if(asymmetric==F){asym <- "symm"}
	if(sampling=="SamplePost_SampleHist"){samp<-"samp_"}
	if(sampling=="FullPost_SampleHist"){  samp<-"fpsh_"}
	if(sampling=="FullPost_FullHist"){    samp<-"fpfh_"}
	if(search=="GPUCBPE"){srch<-"gpucb_"}
	if(search=="Random"){ srch<-"rando_"}
	if(search=="Grid"){   srch<-"grids_"}
	meth <- "elgpriors"; if(method=="Uniform"){meth<-"unipriors"}
	
	info <- sprintf("%s%s%s_%02i_%05s_%s_", samp, srch, sobol, n_init, n_samples, meth)
	info <- paste0(info,paste(models_used, collapse = ''),"_",top_model,"_", seed)
	out <- paste0("DATA_OUTPUT/",info,".csv")
	if(dtype=="FIG"){
		out <- paste0("FIGS/",info,".png")	
	}
	return(out)
}


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
    require(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This is does the summary; it's not easy to understand...
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun= function(xx, col, na.rm) {
                           c( N    = length2(xx[,col], na.rm=na.rm),
                              mean = mean   (xx[,col], na.rm=na.rm),
                              sd   = sd     (xx[,col], na.rm=na.rm),
                              sum  = sum    (xx[,col], na.rm=na.rm)
                              )
                          },
                    measurevar,
                    na.rm
             )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean"=measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}