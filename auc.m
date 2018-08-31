function area=auc(x,y)

if size(x)~=size(y)
    area=NaN;
    warning('The input vectors should have the same size');
    return
else
    n=length(x);
end

[x_, xix]=sort(x,'ascend');
y_=zeros(1,n);
for i=1:n
    y_(i)=y(xix(i));
end
% x__=x_;
% y__=y_;
% for i=1:n
%     if isnan(y_(i)) || isnan(x_(i))
%         y__(i)=[];
%         x__(i)=[];
%     end
% end
x_(isnan(y_))=[];
y_(isnan(y_))=[];
x_(isnan(x_))=[];
y_(isnan(x_))=[];
area=trapz(x_,y_);

