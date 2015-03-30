function fh = logtail_plot(t,logp,logp_cheb,logp_bern,logp_norm,show_legend)
	if nargin < 6
		show_legend = false;
	end
	fh = figure; 
	hold on;
	
	set(0,'DefaultAxesColorOrder',[0 0 0]);
	set(0,'DefaultLineLineWidth',1.5);
	if show_legend
		linewidth = 2;
		plot(nan,'-','linewidth',linewidth);
		plot(nan,'-x','MarkerSize',8,'linewidth',linewidth);
		plot(nan,'--','linewidth',linewidth);
		plot(nan,'-.','linewidth',linewidth);
		legend({'Simulation','Chebyshev','Bernstein','CLT'},...
			'location','south','fontsize',18);
	end
	plot(t,logp,'-');
	plot(t,logp_cheb,'-');
	plot(t(1:5:end),logp_cheb(1:5:end),'x','MarkerSize',5);
	plot(t,logp_bern,'--');
	plot(t,logp_norm,'-.');
	xlim([-abs(min(t))*1.05,max(t)*1.05])
	ylim([min(logp(~isinf(logp)))*1.2,0]);
	xlabel('$$t$$','interpreter','latex','FontSize',18);
	ylabel('$$\widehat{L}(t)$$','interpreter','latex','FontSize',18);
	set(gca,'FontSize',12)
	hold off;
	
end