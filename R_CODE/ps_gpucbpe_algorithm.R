##############################################
###                                        ###
###   Adaptive Design of Experiments for   ###
###       Maximizing Information Gain      ###
###        code by: Brennan Klein          ###
###                                        ###
##############################################
# Parameter-sampling

### Step 1: Sample parameters according to 'uniform' or 'posteriors' from data ### 
get_params <- function(PI, model_num, posteriors, uniform=T){
	# Get one list of three parameters for the specific simulation.
	e <- seq(0,.99,length.out=34)
	a <- seq(0,.99,length.out=34)
	d <- seq(0,.18,length.out=7)
	if(uniform==T){
		epsilon <- sample(e, 1)
		alpha   <- sample(a, 1)
		delta   <- sample(d, 1)
		pi_per  <- sample(c(-1,1),1)*delta + PI
		# epsilon <- sample(e, 1)
		# alpha   <- sample(a, 1)
		# delta   <- sample(d, 1)
		# pi_per  <- sample(c(-1,1),1)*delta + PI
		return(list(epsilon, alpha, pi_per))
	}
	else { # just in case the uniform sampling isn't appropriate 
		a_post <- posteriors[posteriors$variable == "alpha"   & posteriors$model == model_num,]$posterior
		d_post <- posteriors[posteriors$variable == "delta"   & posteriors$model == model_num,]$posterior
		e_post <- posteriors[posteriors$variable == "epsilon" & posteriors$model == model_num,]$posterior
		epsilon <- sample(e, 1, prob=e_post)
		alpha   <- sample(a, 1, prob=a_post)
		delta   <- sample(d, 1, prob=d_post)
		pi_per  <- sample(c(-1,1),1)*delta + PI
		return(list(epsilon, alpha, pi_per))
	}
}

simOneGame <- function(model, A, PI, model_num, sampling_method, posteriors, n_rounds, n_pairs) {
	# simulates a dataset collected from one game
	uniform<-F
	if(sampling_method=="Uniform"){uniform<-T}
	params <- get_params(PI, model_num, posteriors, uniform)
	epsilon <- params[[1]]
	alpha   <- params[[2]]
	pi_per  <- params[[3]]
	
	p1data <- ""
	for(t in c(1:n_rounds)){
		et <- epsilon * exp(- alpha * t)
		pairing <- as.character(sample(seq(1,n_pairs), 1, prob=rep(1/n_pairs,n_pairs)))
		world <- sample(c("a", "b"), 1, prob=c(PI, 1-PI)) 
		# get probability of G or S (depending on playerA or B) in world = world
		p_playerA <- model(pi_per, A, et, alpha, world, p1data)
		p_playerB <- model(pi_per, A, et, alpha, "q", p1data)
		# then, sample actual action
		playerA <- sample(c("G", "S"), 1, prob=c(p_playerA, 1-p_playerA))
		playerB <- sample(c("L", "R"), 1, prob=c(p_playerB, 1-p_playerB))
		game11 <- paste(world, playerA, playerB, sep="")

		if(t == n_rounds){
			p1data <- paste(p1data, game11, sep="")			
		} else{
			p1data <- paste(p1data, game11, "-", sep="")
		}
	}
	dataset <- p1data
	return(dataset)
}

### Step 2: Simulate n games ### 
simMultipleGames <- function(n, model, A, PI, model_num, sampling_method, posteriors, n_rounds, n_players) {
	# Simulates one game, and runs this n times
	if(n_players%%2 !=0){
		print("HEY YOU NEED n_players TO BE EVEN")
	}
	n_pairs <- n_players/2
	
	dat <- replicate(n, simOneGame(model, A, PI, model_num, sampling_method, posteriors, n_rounds, n_pairs))
	return(dat)
}

### Step 3: Take the divergence between the (normalized) history histograms ###
getInfoNumber <- function(A, PI, sampling_method, posteriors, n_samples, asymmetric, sampling, models_used, n_players, n_rounds, current_top_model=1) {
	if(sampling == "SamplePost_SampleHist"){
		allObservedGames <- c()
		mydata2 <- data.frame()
		for(model_num in models_used){
			if(model_num==1){
				m1 <- simMultipleGames(n_samples, model_1n, A, PI, model_num, sampling_method, posteriors, n_rounds, n_players)
				allObservedGames <- unique(c(allObservedGames, m1))
				m1Hist <- data.frame("History"=as.character(m1), "Model"="M1")
				m1Hist$Model <- "M1"
				mydata2 <- rbind(mydata2, m1Hist)
			}
			if(model_num==2){
				m2 <- simMultipleGames(n_samples, model_2n, A, PI, model_num, sampling_method, posteriors, n_rounds, n_players)				
				allObservedGames <- unique(c(allObservedGames, m2))
				m2Hist <- data.frame("History"=as.character(m2), "Model"="M2")
				m2Hist$Model <- "M2"
				mydata2 <- rbind(mydata2, m2Hist)
			}
			if(model_num==3){
				m3 <- simMultipleGames(n_samples, model_3n, A, PI, model_num, sampling_method, posteriors, n_rounds, n_players)				
				allObservedGames <- unique(c(allObservedGames, m3))
				m3Hist <- data.frame("History"=as.character(m3), "Model"="M3")
				m3Hist$Model <- "M3"
				mydata2 <- rbind(mydata2, m3Hist)
			}
			if(model_num==4){
				m4 <- simMultipleGames(n_samples, model_4n, A, PI, model_num, sampling_method, posteriors, n_rounds, n_players)				
				allObservedGames <- unique(c(allObservedGames, m4))
				m4Hist <- data.frame("History"=as.character(m4), "Model"="M4")
				m4Hist$Model <- "M4"
				mydata2 <- rbind(mydata2, m4Hist)
			}
		}
		mydata <- data.frame()
		
		for(model_name in unique(mydata2$Model)){
			model_data <- mydata2[mydata2$Model== model_name,]
			m <- as.character(model_data$History)
			mdf <- data.frame(table( c(m, allObservedGames) )-1)
			colnames(mdf) <- c("History", "Freq")
			mdf$lik <- (mdf$Freq + exp(-16)) / sum(mdf$Freq)
			mdf$Model <- model_name
			mydata <- rbind(mydata, mdf)
		}		
	}
	if(asymmetric==T){ 
		return(asymmetric_divergence(mydata, models_used, current_top_model))
	}
	if(asymmetric==F){
		return(symmetric_divergence(mydata, models_used))
	}
}

##############################################
##############################################

asymmetric_divergence <- function(likelihood_data, models_used, current_top_model=1){
	# Calculates the asymmetric divergence between the models, assuming Model 1 is leading contendor
	models_names <- c("M1", "M2", "M3", "M4")
	model_list <- unique(likelihood_data$Model)
	models <- models_used
	priors <- rep(1/length(model_list), 4)
	prior1 <- priors[as.integer(current_top_model)]
	model_name <- models_names[as.integer(current_top_model)]
	priors_other <- priors[models[models!=current_top_model]]
	non_current_top_models <- model_list[model_list!=model_name]
	hists <- as.character(unique(likelihood_data$History))
	I_sum <- 0
	for(dataset in hists){
		current_compare <- likelihood_data[likelihood_data$History==dataset,]
		l1 <- current_compare[current_compare$Model==model_name,]$lik
		l_others <- current_compare[current_compare$Model!=model_name,]$lik
		denom <- 0 
		for(i in c(1:length(priors_other))){
			denom <- denom + priors_other[i] * l_others[i]
		}
		tmp <- ( (1 - prior1) * l1 ) / denom
	    I_sum <- I_sum + 
	    l1*log(tmp)
	}
	return(I_sum)
}

symmetric_divergence <- function(likelihood_data, models_used=c(1,2,3,4)){	
	Isum <- 0
	for(model_num in models_used){
		Isum <- Isum + asymmetric_divergence(likelihood_data, models_used, model_num)
	}
	I <- Isum/length(models_used)
	return(I)
}
