clc;
clear all;
close all;
file1 = fopen('result.txt','a+');

             
     for load_file = 1:2
    %% initializing variables
    no_part = 10.;
    %% to load file
    switch load_file


       case 1
            file = 'ecoli-0-1_vs_2-3-5';
            test_start =121;
            cvs1=10;
            mus=32;
			epsv=0.4;
       case 2
            file = 'ecoli-0-1_vs_5';
            test_start =121;
            cvs1=10;
            mus=32;
			epsv=0.4;
            
      
        otherwise
            continue;
    end
    % %parameters
         cvs1=[10^-5,10^-4,10^-3,10^-2,10^-1,10^0,10^1,10^2,10^3,10^4,10^5];

         mus=[2^-5,2^-4,2^-3,2^-2,2^-1,2^0,2^1,2^2,2^3,2^4,2^5];

     epsv=[0.1,0.2,0.3,0.4,0.5,0.6];

%Data file call from folder   
filename = strcat('newd/',file,'.txt');
    A = load(filename);
    [m,n] = size(A);
%define the class level +1 or -1    
    for i=1:m
        if A(i,n)==0
            A(i,n)=-1;
        end
    end
% Dividing the data in training and testing    
     test = A(test_start:m,:);
    train = A(1:test_start-1,:);


    [no_input,no_col] = size(train);
  
    x1 = train(:,1:no_col-1);
    y1 = train(:,no_col);
	    
    [no_test,no_col] = size(test);
    xtest0 = test(:,1:no_col-1);
    ytest0 = test(:,no_col);

    A=[x1 y1];    %training data
    A_test=[xtest0,ytest0];    %testing data
    [m,n] = size(A);
    p=0;
    for i=1:m
        if A(i,n)==1
            p=p+1;
        end
    end
    d=(m-2*p);


    [lengthA,n] = size(A);
    min_err = -10^-10.;

     [no_input,no_col]=size(A);
  obs = A(:,no_col);   
    C=A;
    A = [];
 B = [];
 u=abs(d);
for i = 1:no_input
    if(obs(i) == 1)
        A = [A;C(i,1:no_col-1)];
    else
        B = [B;C(i,1:no_col-1)];
    end;
end;

sb1=size(A,1);
sb=size(B,1);
ptb1=sb1/u;
ptb=sb/u;
Au=A(1:ptb1:sb1,:);
Bu=B(1:ptb:sb,:);
di=size(Au,1)-size(Bu,1);
if(di>0)
Bu=[Bu ;Bu(1:abs(di),:)];
elseif(di<0)
Au=[Au ;Au(1:abs(di),:)];
end   
 U=(Au+Bu)/2;   
   A=C;


  for C1 = 1:length(cvs1)
            c = cvs1(C1)

            for mui=1:length(mus)
                mu=mus(mui);

             for  epsi = 1:length(epsv)
                    e = epsv(epsi);
                    avgerror = 0;
                    block_size = lengthA/(no_part*1.0);
                    part = 0;
                    t_1 = 0;
                    t_2 = 0;
                    while ceil((part+1) * block_size) <= lengthA
                   %% seprating testing and training datapoints for
                   % crossvalidation
                                t_1 = ceil(part*block_size);
                                t_2 = ceil((part+1)*block_size);
                                B_t = [A(t_1+1 :t_2,:)];
                                Data = [A(1:t_1,:); A(t_2+1:lengthA,:)];
                   %% testing and training
                                [accuracy_with_zero,time] = rutsvm(Data,B_t,U,c,mu,e);
                                avgerror = avgerror + accuracy_with_zero;
                                part = part+1
                     end
       
           if avgerror > min_err
               min_err = avgerror;
               min_mu = mu;
               min_c1 = c;
                min_e=e;
           end

             end
            end
  end
       
   
 [accuracy,time] = rutsvm(A,A_test,U,min_c1,min_mu,min_e);
 fprintf(file1,'%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\tu=%g\n',file,size(A,1),size(A_test,1),accuracy,min_c1,min_mu,min_e,time,u);
 
end