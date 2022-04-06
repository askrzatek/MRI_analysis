function [result log a]= do_statistique_test(group1,group2)


% normal distribution test shapiro wilk
%alpha =0.05

a=0;
[Hsw1, pValue1, W1] = swtest(group1)
[Hsw2, pValue1, W2] = swtest(group2)

if Hsw1 == 1 || Hsw2==1
  a=1;  
    log = sprintf('   Not normally distributed data: H1= %d, H2= %d \n',Hsw1,Hsw2);
    log = sprintf('%s   Using non-parametric test : Wilcoxon rank sum ''equivalent to the Mann-Whitney'' \n',log)
    
    [p,h,stats] = ranksum(group1,group2)
   
    log = sprintf('%s   Result << H= %d, P= %d >>\n',log,h,p);
    
    result.h =h;
    result.p =p;
    
elseif Hsw1 ==0 & Hsw2 == 0
    a=0;
    
    log = sprintf('   Normally distributed data: H1= %d, H2= %d \n',Hsw1,Hsw2);
    log = sprintf('%s   Using parametric test: comparison of the two variances with Fisher test\n',log)
    
    [hf,pf,cif,statsf] = vartest2(group1,group2)
    
    
    if hf == 0
        log = sprintf('%s   var1 = var2 : <homoscedasticity> so used ttest2\n',log)
        [H,P,CI,STATS] = ttest2(group1,group2)
        log =sprintf('%s   Result : H = %d, P = %d\n',log,H,P)
        result.h=H;
        result.p=P;
        result.ci=CI;
        
        result.stats=STATS;
    else 
        log = sprintf('%s   var1 =/= var2 : <heteroscedasticity> so used Aspin Welich test\n',log)
        [H,P,CI,STATS] = ttest2(group1,group2, 'Vartype','unequal')
        
        log =sprintf('%s   Result : H = %d, P = %d\n',log,H,P)
        
        result.h=H;
        result.p=P;
        result.ci=CI;
        result.stats=STATS;
    end
end
end
    