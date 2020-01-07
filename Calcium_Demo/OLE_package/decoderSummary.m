% Visually Assess Decoding Results...
idx=1:size(f,1);
figure(1)
clf
subplot(2,1,1)
h = fspecial('disk',2);
imagesc(1:size(f,1),pvec,imfilter(f,h)') % imfilter--- multidimensional filtering 

xlim(xl)
hold on
plot(1:length(idx),rpos(idx)*pi,'r.')
[tmp,maxpost]=max(f');
plot(1:length(idx),maxpost/length(pvec)*2*pi,'k.')
hold off
ylabel('Position','fontsize',20)
colormap(flipud(cbrew(256,'rdbu')))
set(gca,'fontsize',20); 

subplot(2,1,2)
plot(mean(spk(idx,:),2),'k')
ylabel('Avg Spike count (MUA)','fontsize',20)
xlabel('Time [bins]','fontsize',20)
xlim(xl)
box off
  set(gca,'fontsize',20);
    

% Average prediction matrix...
edges = [pvec-mean(diff(pvec))/2 max(pvec)+mean(diff(pvec))/2];
[tmp,rpos_bin] = histc(rpos(idx)*pi,edges);
C=[];
for i=1:length(pvec)
    bidx = find(rpos_bin==i);
    C(:,i) = mean(f(bidx,:),1);
end
if all(C(:)>=0)
    C = bsxfun(@rdivide,C,sum(C));
end
figure(2)
clf
imagesc((C))
axis image
colorbar
colormap(flipud(cbrew(256,'rdbu')))
title('Average Linear Decoding')
xlabel('True Position')
ylabel('Predicted Position')
  set(gca,'fontsize',20);
    

% Circular error...
err =  circ_dist(     rpos(idx)*pi, maxpost/length(pvec)*2*pi)*trackLen/pi;
errp = circ_dist(     rpos(idx)*pi, maxpost/length(pvec)*2*pi)*trackLen/pi;
errm = circ_dist(2*pi-rpos(idx)*pi, maxpost/length(pvec)*2*pi)*trackLen/pi;


