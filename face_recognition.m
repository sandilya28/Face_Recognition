clear all;
close all;
names = {'centerlight','glasses','happy','leftlight','noglasses','normal','rightlight','sad','sleepy','surprised','wink'};
face = 0;
set1 = zeros(33,2);
set2=zeros(33,2);
set3=zeros(33,2);
set4=zeros(33,2);
set5=zeros(33,2);
count1 = 1;
for i = 1:165
    if(mod(i,11) == 1)
        a = randperm(11);
        count = 1;
        face = face + 1;
    end 
    if(mod(i,5) == 1)
            set2(count1,1) = face;
            set2(count1,2) = a(count);
            count = count+1;
    elseif(mod(i,5)==2)
            set3(count1,1) = face;
            set3(count1,2) = a(count);
            count = count+1;
    elseif(mod(i,5)==3)
            set4(count1,1) = face;
            set4(count1,2) = a(count);
            count = count+1;
    elseif(mod(i,5)==4)
            set5(count1,1) = face;
            set5(count1,2) = a(count);
            count = count+1;
    elseif(mod(i,5)==0)
            set1(count1,1) = face;
            set1(count1,2) = a(count);
            count = count+1;
            count1 = count1 + 1;
    end
end        

Data = [];
for i = 1:5
% finding the total Data matrix from the 5 sets
    if(i==1)
        b = set1;
    elseif(i==2)
        b = set2;
    elseif(i==3)
        b = set3;
    elseif(i==4)
        b = set4;
    elseif(i==5)
        b = set5;
    end
    for j = 1:33
% for each row in a set
        if b(j,1) <= 9
            filename = strcat('subject','0',num2str(b(j,1)),'.',names(b(j,2)));        
        else
            filename = strcat('subject',num2str(b(j,1)),'.',names(b(j,2)));
        end
        
        filename = sprintf('./yalefaces/%s',filename{:});
        matrix = double(imread(filename));
        row = transpose(matrix(:));
        Data = [Data;row];
    end
end

% 5-fold verification
for i=1:5
    % Test data and Training Data
    Data1 = [];
    Test = [];
    if i==1 
        for j = 1:33
            Test = [Test;Data(j,:)];
        end
        for j = 34:165
            Data1 = [Data1;Data(j,:)];
        end
    elseif i==2
        for j = 1:33
            Data1 = [Data1;Data(j,:)];
        end
        for j = 34:66
            Test = [Test;Data(j,:)];
        end
        for j = 67:165
            Data1 = [Data1;Data(j,:)];
        end
    elseif i==3
        for j = 1:66
            Data1 = [Data1;Data(j,:)];
        end
        for j = 67:99
            Test = [Test;Data(j,:)];
        end
        for j = 100:165
            Data1 = [Data1;Data(j,:)];
        end
    elseif i==4
        for j = 1:99
            Data1 = [Data1;Data(j,:)];
        end
        for j = 100:132
            Test = [Test;Data(j,:)];
        end
        for j = 133:165
            Data1 = [Data1;Data(j,:)];
        end
    elseif i==5
        for j = 1:132
            Data1 = [Data1;Data(j,:)];
        end
        for j = 133:165
            Test = [Test;Data(j,:)];
        end
    end
    
    % Average over the features
    Average = sum(Data1);
    Average=Average/132;
    xxx = size(Average);
    T = ones(132,xxx(2));
    % Calculating Mean centered Data Matrix T
    for j = 1:132
        for k = 1 : xxx(2)
            T(j,k) = Data1(j,k) - Average(k);
        end
    end
    
    % S is the covariance Matrix
    S = T * transpose(T);
    
    % Calculating eigen vectors for S
    [e_vec,e_val]=eig(S);
%     e_vec = transpose(T)*e_vec;
%     e_vec = transpose(e_vec);
%     e_val = transpose(T)*e_val;
%     e_val = transpose(e_val)
    
    e_val = transpose(diag(e_val));
    e_val = fliplr(e_val);
    
    % Threshold 0.95 
    e_sum = sum(e_val);
    temp_sum = 0.0;
    
    threshold = 1;
    while(temp_sum/e_sum < 0.95)
        temp_sum = temp_sum + e_val(threshold);
        threshold = threshold + 1;
    end
    threshold = threshold - 1;
    
    %removing unnecessary eigen vectors
    e_vec = fliplr(e_vec);
    for j = threshold+1:132
        e_vec(:,threshold+1)=[];
    end
     
    e_vec = transpose(T)*e_vec;
    e_vec = transpose(e_vec);
    
    si=[];
    for j=1:132
        si = [si ; transpose(double(e_vec) * double(transpose(T(j,:))))];
    end
    
end

        
        
        