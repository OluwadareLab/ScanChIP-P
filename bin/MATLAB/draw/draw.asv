
Chr_Data = dlmread('sym_30x30.txt'); 

Visualize(Chr_Data);

Tadlist = dlmread('tad.txt');

for i = 1:length(Tadlist(:,1))
     hold on;
     Start = Tadlist(i,1);
     Last = Tadlist(i,2);
     for j = Start:Last
             plot(Start,j,'r.','MarkerSize',10);
             plot(j,Start,'r.','MarkerSize',10);
             plot(j,Last,'r.','MarkerSize',10);
             plot(Last,j,'r.','MarkerSize',10);
     end        
end
fprintf('The Number of TAD = %d\n', length(Tadlist(:,1)));
