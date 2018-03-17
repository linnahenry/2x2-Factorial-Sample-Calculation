




#this is for a 2^2 factorial 
 
effect_sizes <- c(5, 6, 10, 11) #for treatments: aa, ab, ba, bb respectively
effect_var <- 5 #if pilot study variance is unknown, function can be run multiple times with different variances
potential_n <- seq(20, 400, by = 20) #n per group -- for total n you multiply n*4
alpha_level <- 0.05

####Run this more than once (simulations) & take average. 

calc_power <- function(effect_sizes, effect_var, potential_n, alpha_level = 0.05){
	n_out <- c()
	pa_out <- c()
	pb_out <- c()
	pab_out <- c()
	fa_out <- c()
	fb_out <- c()
	fab_out <- c()
	for(i in 1:length(potential_n)){
		dq11 <- rnorm(potential_n[i], mean = effect_sizes[1], sd = sqrt(effect_var))
		dq12 <- rnorm(potential_n[i], mean = effect_sizes[2], sd = sqrt(effect_var))
		dq21 <- rnorm(potential_n[i], mean = effect_sizes[3], sd = sqrt(effect_var))
		dq22 <- rnorm(potential_n[i], mean = effect_sizes[4], sd = sqrt(effect_var))
		y11 <- sum(dq11)
		y12 <- sum(dq12)
		y21 <- sum(dq21)
		y22 <- sum(dq22)
		
		yijk <- y11 + y12 + y21 + y22
		nn <- potential_n[i]
		n_out <- c(n_out, nn*4)
		ssa <- (((y12+y11)**2 + (y22+y21)**2)/(2*(nn))) - (((yijk)**2)/(4*(nn)))
		ssb <- (((y12+y22)**2 + (y11+y21)**2)/(2*(nn))) - (((yijk)**2)/(4*(nn)))
		
		yij2 <- sum(dq11**2) + sum(dq12**2) + sum(dq21**2) + sum(dq22**2)
		
		sssub <- ((y11**2 + y12**2 + y21**2 + y22**2)/nn) - ((yijk**2)/(nn*4))
		
		ssab <- sssub - ssa - ssb
		
		sst <- yij2 - ((yijk**2)/(nn*4))
		
		sse <- sst - sssub
		
		msa <- ssa
		msb <- ssb
		msab <- ssab
		mse <- sse/(4*((nn)-1))
		
		f_a <- msa/mse
		f_b <- msb/mse
		f_ab <- msab/mse
		
		fa_out <- c(fa_out, f_a)
		fb_out <- c(fb_out, f_b)
		fab_out <- c(fab_out, f_ab)
		
		fcrit <- qf(alpha_level, 1, nn-1, lower.tail = FALSE)
		
		f_a_power <- 1-pf(fcrit, 1, nn-1, f_a)
		f_b_power <- 1-pf(fcrit, 1, nn-1, f_b)
		f_ab_power <- 1-pf(fcrit, 1, nn-1, f_ab)
		
		pa_out <- c(pa_out, f_a_power)
		pb_out <- c(pb_out, f_b_power)
		pab_out <- c(pab_out, f_ab_power)
		
	}
	output <- cbind(n_out, fa_out, pa_out, fb_out, pb_out, fab_out, pab_out)
	output <- data.frame(output)
	colnames(output) <- c("n total", "F-val Factor A", "Power Factor A", "F-val Factor B", "Power Factor B", "F-val Factor A*B", "Power Factor A*B")
	return(output)		
}	

calc_power(effect_sizes, effect_var, potential_n, alpha_level = 0.05)




