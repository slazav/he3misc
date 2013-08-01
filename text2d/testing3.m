load('/rota/Analysis/NMRcalc/spinwaves/TextLib_Slava1/testing_gr.dat')
load('/rota/Analysis/NMRcalc/spinwaves/TextLib_Slava/testing_gr_original.dat')

testing_gr(18:length(testing_gr),:)=[];
testing_gr_original=testing_gr_original(1:5,:);

%check if gradient is the same 
inconsistency_perc=abs((sum(testing_gr_original(:,6))-10/(2*pi)*sum(testing_gr(:,11)))/sum(testing_gr_original(:,6)))

%sum(testing_gr([1,2,5,6],5).*testing_gr([1,2,5,6],6))/2

sum(testing_gr([1,2,5,6],7))/2
testing_gr_original(2,4)


new=sum(testing_gr(:,6)+testing_gr(:,9));