function [axes1,figure1]=figureX1(size)
if nargin<1
    size=[7,4.33];
end
figure1=figure('InvertHardcopy','off','Units','centimeters','Position',[25,10,size],'Color',[1 1 1],'Renderer','painters');
axes1 = axes('Parent',figure1);
% set(axes1,'FontName','times new roman','FontSize',10);

end