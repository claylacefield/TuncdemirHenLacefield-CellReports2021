function [RateCorr]= rateCorrelation (Matched, matrix)
%calculates rate correlation between position vs rate matrix from two sessions 
%that are 'Matched' across sessions 
% 'matrix' is session numbers
popVecCorr = corr(Matched{matrix(1)}, Matched{matrix(2)}, 'rows', 'pairwise' );
cellCorr = corr(Matched{matrix(1)}', Matched{matrix(2)}', 'rows', 'pairwise'); 

ShuffpopVecCorr=[];
ShuffcellCorr=[];
   rng('shuffle')
   for sh = 1:200
       posRatesCellRand = Matched;
       posRatesCellRand{matrix(1)} = Matched{matrix(1)}(randperm(size(Matched{matrix(1)}, 1)), :);
       posRatesCellRand{matrix(2)} = Matched{matrix(2)}(randperm(size(Matched{matrix(2)}, 1)), :);
ShuffpopVecCorr{sh} = corr(posRatesCellRand{matrix(1)}', Matched{matrix(2)}', 'rows', 'pairwise' );
ShuffcellCorr{sh}= corr(posRatesCellRand{matrix(1)}, Matched{matrix(2)}, 'rows', 'pairwise' );
   end
    ShpopVecCorr=[]; 
   ShcellCorr=[];
   for i=1:length(ShuffpopVecCorr)
        a=ShuffpopVecCorr{i};
       c=ShuffcellCorr{i};
       a=diag(a); 
       c=diag(c); 
        ShpopVecCorr=[ShpopVecCorr a];
       ShcellCorr=[ShcellCorr c];
   end
   
  RateCorr.popVecCorr=popVecCorr;
  RateCorr.cellCorr=cellCorr;
  RateCorr.ShpopVecCorr=ShpopVecCorr;
  RateCorr.ShcellCorr=ShcellCorr;
