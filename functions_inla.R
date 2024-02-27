################################################################################
# Useful functions to facilitate analysis
# - Created on 13-01-2021
################################################################################

## Get the fixed effect estimates as exp(beta) (for binom = odds)
get_fixed_effects_inla <- function(fit_out,OR=TRUE){
  
  # Get fixed on logit-level
  tidy.inla <- function(x){
    
    # x = model_inla
    term_names <- rownames(x$summary.fixed)
    
    tibble::as_tibble(x$summary.fixed) %>%
      dplyr::mutate(terms = term_names) %>%
      dplyr::select(param = terms,
                    dplyr::everything())
  }
  
  if(!OR){
    
    marginals <- tidy.inla(fit_out) %>% 
      rename(lower = `0.025quant`,
             upper = `0.975quant`,
             median = `0.5quant`) %>% 
      select(-c("median","mode","kld"))
    
    
  } else{
    
    marginals <- data.frame(matrix(NA,nrow=length(fit_out$marginals.fixed),ncol=5))
    marginals[,1] <-names(fit_out$marginals.fixed)
    
    for (i in 1:length(fit_out$marginals.fixed)){
      
      m_or <-inla.tmarginal(fun=function(x) exp(x), fit_out$marginals.fixed[[i]])
      
      est <- inla.zmarginal(m_or,silent=TRUE)
      marginals[i,2] <- est$mean
      marginals[i,3] <- est$sd
      
      # summary marginal distribution: HPD interval
      hpd <- inla.hpdmarginal(.95,marginal=m_or)
      marginals[i,4] <- hpd[,1]
      marginals[i,5] <- hpd[,2]
      
    }
    names(marginals) = c("param","mean","sd","lower","upper")
  }
  
  # for (i in 1:length(fit_out$marginals.fixed)){
  #   
  #   if(OR){
  #     m_or <-inla.tmarginal(fun=function(x) exp(x), fit_out$marginals.fixed[[i]])
  #     
  #   }else{
  #     m_or <-fit_out$marginals.fixed[[i]]
  #     
  #   }
  #   
  
  return(marginals)
}



# Get Random effects (exp(effect))
get_random_effects_inla <- function(inla.model,OR=TRUE,nis_ids=NULL){
  
  # Binomial model -> odd ratios
  
  # empty list 
  re_estimates <- list()
  
  for(re in names(inla.model$marginals.random)){
    
    
    marginals_re <- data.frame(matrix(NA,nrow=length(inla.model$marginals.random[[re]]),ncol=5))
    marginals_re[,1] <-c(names(inla.model$marginals.random[[re]]))
    
    for (i in 1:length(inla.model$marginals.random[[re]])){
      
      if(OR){
        marg <-inla.tmarginal(fun=function(x) exp(x), inla.model$marginals.random[[re]][[i]])
        
      }else{
        marg <-inla.model$marginals.random[[re]][[i]]
        
      }
      
      est <- inla.zmarginal(marg,silent=T)
      marginals_re[i,2] <- est$mean
      marginals_re[i,3] <- est$sd
      
      # 95% CI (HPD interval)
      hpd <- inla.hpdmarginal(.95,marginal=marg)
      marginals_re[i,4] <- hpd[,1]
      marginals_re[i,5] <- hpd[,2]
      
    }
    
    # Set col names
    names(marginals_re) = c("param","mean","sd","lower","upper")
    
    if(re=="practice_id"){
      
      marginals_re <- marginals_re %>% 
        mutate(param=inla.model$summary.random$practice_id$ID)
      
      
    } else {
      
      if(is.null(nis_ids)){
        
        marginals_re <- marginals_re %>%
          mutate(param=inla.model$summary.random[[re]]$ID)
        
      } else 
        
        marginals_re <- marginals_re %>%
          mutate(param=inla.model$summary.random[[re]]$ID) %>% 
          left_join(nis_ids,by="ID") %>% 
          mutate(province = case_when(
            substr(nis_code,1,1)==1 ~ "Antwerpen",
            substr(nis_code,1,1)==2 ~ "Vlaams-Brabant",
            substr(nis_code,1,1)==3 ~ "West-Vlaanderen",
            substr(nis_code,1,1)==4 ~ "Oost-Vlaanderen",
            TRUE ~ "Limburg")
          ) %>% 
          mutate(province=factor(province,levels=c("West-Vlaanderen","Oost-Vlaanderen",
                                                   "Vlaams-Brabant","Antwerpen",
                                                   "Limburg")))
    }
    
    # Append to re_estimates list
    re_estimates <- append(re_estimates,list(marginals_re))
  }
  
  
  names(re_estimates)=names(inla.model$marginals.random)
  return(re_estimates)
  
}

#-----------------#-----------------#-----------------#-----------------#-------

# Get Goodness of Fit statistics
get_model_fit_statistics_inla <- function(fit_out){
  
  marginals <- data.frame(matrix(NA,nrow=1,ncol=5))
  names(marginals) <-c("WAIC","DIC","CPO","MMLIK","Eff no par (WAIC)")
  # names(marginals) <-c("WAIC","DIC","CPO","Eff no par (WAIC)")
  
  marginals[1,1] <- fit_out$waic$waic
  marginals[1,2] <- fit_out$dic$dic
  marginals[1,3] <- sum(log(fit_out$cpo$cpo),na.rm=T)
  marginals[1,4] <- fit_out$misc$configs$max.log.posterior
  marginals[1,5] <- fit_out$waic$p.eff
  
  return(marginals)
  
}

#-----------------#-----------------#-----------------#-----------------#-------

# Function to plot fixed effects 
plot_fixed_effects <- function(fe_estimates,OR=TRUE){
  
  if(OR){
    x_intercept=1
  }else{
    x_intercept=0
  }
  
  if("model"%in%colnames(fe_estimates)){
    
    ggplot(data = fe_estimates,
           mapping = aes(y=param,x=mean,color=model,shape=model,linetype=model))+ 
      geom_point()+
      geom_errorbar(aes(xmin = lower, xmax = upper),width=0,size=1)+
      geom_vline(xintercept=x_intercept)+
      labs(x="Odds ratio",title="Fixed effects estimates")+
      scale_linetype_manual(values=c("solid","dashed","11"))+
      scale_x_continuous(breaks=seq(round(min(fe_estimates$mean),0),round(max(fe_estimates$mean),0),0.5))+
      theme(legend.position = "bottom")
    
  } else {
    
    ggplot(data = fe_estimates %>% 
             mutate(param = factor(case_when(
               param == "(Intercept)" ~ "Intercept",
               param == "age_category(-1,5]" ~ "Age group: 0-5",
               param == "age_category(5,16]" ~ "Age group: 6-16",
               param == "age_category(50,65]" ~ "Age group: 51-65",
               param == "age_category(65,85]" ~ "Age group: 66-85",
               param == "age_category(85,110]" ~ "Age group: 85+",
               param == "increased_compensationYes" ~ "Increased compensation",
               param == "sexM" ~ "Male"
             ), 
             levels = c("Intercept","Age group: 0-5","Age group: 6-16",
                        "Age group: 51-65","Age group: 66-85","Age group: 85+",
                        "Increased compensation","Male")
             )
             ),
           mapping = aes(y = param,x = mean))+ 
      geom_point()+
      geom_errorbar(aes(xmin = lower, xmax = upper),width=0,size=1)+
      geom_vline(xintercept = x_intercept)+
      labs(
        x="Posterior mean estimate",
        y=""
        #title="Fixed effects estimates"
      )+
      scale_linetype_manual(values=c("solid","dashed","11"))+
      theme_bw()+
      scale_x_continuous(breaks=seq(round(min(fe_estimates$mean),0),round(max(fe_estimates$mean),0),0.5))+
      theme(legend.position = "bottom",
            panel.grid = element_blank())
    
  }
}

#-----------------#-----------------#-----------------#-----------------#-------

# Function to plot random effect 
plot_random_effects <- function(re_estimates, n_re=2, map=TRUE, plot_map){
  
  ## 31/08: change left join by ID to param (because changed in get_....effects)
  
  if(map){
    
    # # Load map as SF object 
    # source("/Users/u0121893/OneDrive - KU Leuven/PhD/Research/Simulation/Spatial/Code/functions_geo.R")
    # map_sf.shp <- get_sf_map_flanders()
    # map_sf.shp <- map_sf.shp %>% 
    #   arrange(NISCODE)
    # map_sf.shp$ID=1:nrow(map_sf.shp)
    
    map_sf.shp <- plot_map
    
    pre<-list()
    
    for(re in names(re_estimates)[1:n_re]){
      
      if(re == "practice_id"){
        
        # load pracitce locations 
        practice_coordinates <- read.csv("/Users/u0121893/OneDrive - KU Leuven/phd/projects/simulation/spatial/github-repo-data2/data/practice_location.csv")
        practice_re <- practice_coordinates %>% 
          left_join(re_estimates$practice_id,by=c("practice_id"="param"))
        
        
        pre[[re]] <- ggplot(data = map_sf.shp) +
          theme_void() +
          geom_sf(fill = 'grey85', color = 'gray', alpha = .1) +
          geom_point(data = practice_re, aes(x = lon.pc_mean, y = lat.pc_mean,size=abs(mean),colour=mean), alpha = 0.8)+
          theme(legend.position = "right",
                plot.margin = margin(.1,.1,.1,.1, unit="cm"),
                text = element_text(size=18)) + 
          scale_colour_continuous(expression(v(z)),type="viridis") +
          # scale_colour_gradient(expression(V(z)),low="grey80",high="grey10") + 
          scale_size(guide='none') 
        # labs(title="Practice Random Effect (Odds ratio)")
        
        
        
      }else if(re=="ID") {
        
        pre[[re]] <- ggplot(data = map_sf.shp %>% left_join(re_estimates[[re]],by=c("ID"="param")))+
          theme_void()+
          geom_sf(aes(fill=mean))+
          scale_fill_continuous(expression(u(x)),type="viridis", n.breaks = 6)+
          # scale_fill_gradient(expression(U(x)),low="grey90",high="grey10") + 
          theme(plot.margin = margin(.1,.1,.1,.1, unit="cm"),
                text = element_text(size=18))
        #labs(title="Community Random Effect (Odds ratio)")
      }else{
        
        pre[[re]] <- ggplot(data = map_sf.shp %>% left_join(re_estimates[[re]],by=c("ID"="param")))+
          theme_void()+
          geom_sf(aes(fill=mean))+
          scale_fill_continuous("Spatial",type="viridis")#+
        #labs(title="Community Random Effect (Odds ratio)")
        
        
      }
      
      # library(gridExtra)
      library(cowplot)
      plots_re<-plot_grid(plotlist = pre,ncol=1,align="hv",axis="r")
      
    }
    
  } else {
    
    plots <-list()
    
    for(re in names(re_estimates)[1:n_re]){
      
      if(re=="practice_id"){
        
        plots[[re]] <- ggplot(data=re_estimates[[re]] %>% mutate(ID = as.numeric(factor(ID)))) +
          geom_line(aes(x=ID,y=mean))+
          geom_ribbon(aes(x=ID,ymin=lower,ymax=upper),alpha=0.2,linetype=2)+
          labs(title=paste0("Random Effect: ",re))
        
      } else {
        
        if("province" %in% colnames(re_estimates[[re]])){
          
          plots[[re]] <- ggplot(data=re_estimates[[re]] %>% 
                                  arrange(province,nis_code) %>% 
                                  mutate(ID = 1:nrow(re_estimates[[re]]))) +
            geom_point(aes(x=ID,y=mean,colour=province))+
            geom_errorbar(aes(x=ID,ymin=lower,ymax=upper,colour=province),alpha=0.2,linetype=2)+
            labs(title=paste0("Random Effect: ",re))
          
        }else{
          plots[[re]] <- ggplot(data=re_estimates[[re]]) +
            geom_line(aes(x=ID,y=mean))+
            geom_ribbon(aes(x=ID,ymin=lower,ymax=upper),alpha=0.2,linetype=2)+
            labs(title=paste0("Random Effect: ",re))
        }
      }
    }
    
    library(gridExtra)
    plots_re<-do.call(grid.arrange,plots)
    
  }
  
  return(plots_re)
}
#-----------------#-----------------#-----------------#-----------------#-------


## Get Predicted number of cases from binom moel
get_predicted_incidence_inla <- function(inla_model,model_data){
  
  # PREDICTIONS
  predicted.values.mean <- c()
  predicted.values.std <- c()
  for(i in 1:length(inla_model$marginals.fitted.values)){
    predicted.values.mean[i] <- 
      inla.emarginal(function(x) x*model_data$npatients[i],
                     inla_model$marginals.fitted.values[[i]])
    
    # predicted.values.std[i]<-
    #   inla.zmarginal(model_data$npatients[i]*
    #                    inla_model$marginals.fitted.values[[i]],silent=T)[[2]]
    
  }
  
  model_data$fitted.values <- round(predicted.values.mean, digits=3)
  # model_data$fitted.values.std <- round(predicted.values.std, digits=3)
  data_fit <- model_data[is.na(model_data$ncases),]
  
  #RR.mean represents the predicted number of cases in that area
  data_fit <- data_fit %>% 
    group_by(nis_code) %>% 
    summarise(RR.mean=sum(fitted.values),Totpop=sum(npatients))
  
  #incidence in a population of 100,000 inhabitants
  data_fit$pred <- (data_fit$RR.mean/data_fit$Totpop)*100000
  
  return(data_fit)
  
}
#-----------------#-----------------#-----------------#-----------------#-------



# from utils.R
Fisher_bin <- function(x, num_bin = 5) {
  require(classInt)
  class <- classIntervals(x[!is.na(x)], num_bin, style = 'fisher')
  ret <- cut(x, breaks = class$brks, right = FALSE, include.lowest = TRUE, dig.lab = 6)
  return(ret)
}


## Plot predicted incidence
plot_predicted_incidence <- function(predictions_inla,nclr=8){
  
  plot_breaks <- Fisher_bin(predictions_inla$pred,nclr)
  predictions_inla$incidence <- plot_breaks
  plotclr <- rev(brewer.pal(nclr,"RdYlGn"))
  
  source("Z:/Projects_Arne/Simulation_Spatial/Code/functions_geo.R")
  map_sf.shp <- get_sf_map_flanders()
  map_sf.shp <- map_sf.shp %>% 
    arrange(NISCODE)
  
  ggplot(data = map_sf.shp %>% 
           left_join(predictions_inla,by=c("NISCODE"="nis_code")))+
    theme_void()+
    geom_sf(aes(fill=incidence))+
    scale_fill_manual("LRTI cases in 100.000 inhabitants",values=plotclr)
  
}
#-----------------#-----------------#-----------------#-----------------#-------
