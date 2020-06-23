function [accuracy,time]=rutsvm(C,test_data,U,c1,mu,epsilon)
  [no_input,no_col]=size(C);
  obs = C(:,no_col);    
 A = [];
 B = [];

for i = 1:no_input
    if(obs(i) == 1)
        A = [A;C(i,1:no_col-1)];
    else
        B = [B;C(i,1:no_col-1)];
    end
end
d=abs(size(A,1)-size(B,1));
[ux uy]=size(U);

U=U(1:(ux-1)/d:(ux-(ux-1)/d),:);

c2=c1;
c3=c1;
       
    
    ep = 0.00001;

    [m3,n] = size(U);
     e3 = ones(m3,1);
     
     [m1,n] = size(A); 

     cl=ceil(m1/2);
          
    e1 = ones(m1,1); 
    [m2,n] = size(B);
   
     e2 = ones(m2,1);
    
     B1=B(1:(m2-1)/m1:(m2-(m2-1)/m1),:);
    
      ed=ones(d,1);

    C = [A ; B1];  
     m= 2*m1;
     tic

     K=zeros(m1,m);
     for i=1:m1
        for j=1:m
            nom = norm( A(i,:)  - C(j,:)  );
            K(i,j) = exp( -1/(2*mu*mu) * nom * nom );
        end
    end
       
     H = [K e1];
         
    
     K=zeros(m2,m);
    for i=1:m2
        for j=1:m
            nom = norm( B(i,:)  - C(j,:)  );
            K(i,j) = exp( -1/(2*mu*mu) * nom * nom );
        end
    end

    G = [K e2];
    

     K=zeros(m3,m);
     for i=1:m3
        for j=1:m
            nom = norm( U(i,:)  - C(j,:)  );
            K(i,j) = exp( -1/(2*mu*mu) * nom * nom );
        end
    end

    O = [K e3];
    

    em = m+1;
    
  lowb1=zeros(m1+cl,1);
    lowb2=zeros(m1+d,1);
    
   upb1 = [c2*e1(1:m1);c3*ed(1:cl)];
   upb2 = [c2*e1;c3*ed];
   g=(m2-1)/m1;
   g1=m2-g;
   ua=(m3-1)/cl;
   ua1=m3-ua;
     

    HTH = H' * H;
    invHTH = inv(HTH + ep * speye(em) );

    G1=[G(1:g:g1,:);-O(1:ua:ua1,:)];
    GOINVGOT = G1 * invHTH * G1';
    
    GTG = G' * G;
    invGTG = inv (GTG + ep * speye(em));
    HO=[H;O];
    HOINVHOT = HO * invGTG * HO';
 

       f1 = -[e2(1:m1);(epsilon-1)*ed(1:cl)]';
       f2 = -[e1;(1-epsilon)*ed]';
    
    u1 = quadprog(GOINVGOT,f1,[],[],[],[],lowb1,upb1);
    u2 = quadprog(HOINVHOT,f2,[],[],[],[],lowb2,upb2);
    time= toc

    w1= - invHTH * (G1(1:m1,:)'*u1(1:m1)-O(1:ua:ua1,:)' *u1(m1+1:m1+cl));
    w2 =  invGTG * (H' *u2(1:m1)+O'*u2(m1+1:m1+d));
    [no_test,no_col] = size(test_data);   

     for i=1:no_test
        for j=1:m
            nom = norm( test_data(i,1:no_col-1)  - C(j,:)  );
            Ker_row(i,j) = exp( -1/(2*mu*mu) * nom * nom );
            
        end
     end
     K = [Ker_row ones(no_test,1)];
   
     y1 = K * w1 / norm(w1(1:size(w1,1)-1,:));
     y2 = K * w2 / norm(w2(1:size(w2,1)-1,:));
     
    for i = 1 : no_test
    if abs(y1(i)) < abs(y2(i))
        classifier(i) = 1;
    else
        classifier(i) = -1;
    end;
end;
%-----------------------------
match = 0.;
match1=0;
classifier = classifier';
obs1 = test_data(:,no_col);
posval=0;
negval=0;

for i = 1:no_test
    if(obs1(i)==1)
    if(classifier(i) == obs1(i))
        match = match+1;
    end
     posval=posval+1;
    elseif(obs1(i)==-1)
        if(classifier(i) ~= obs1(i))
        match1 = match1+1;
        end
    negval=negval+1;
    end
end
if(posval~=0)
a_pos=(match/posval)
else
a_pos=0;
end

if(negval~=0)
am_neg=(match1/negval)
else
am_neg=0;
end

AUC=(1+a_pos-am_neg)/2;

accuracy=AUC*100;
