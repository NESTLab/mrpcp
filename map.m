%% Visualization of the path generated for each robot

function [] = map(A, robot, T, D, K, N, x_pos, y_pos, L)
    figure;
    G = digraph(A);
    h = plot(G,'XData', x_pos, 'YData', y_pos);
    axis square
    hold on
    title("Path of Robot "+num2str(robot));
    highlight(h,[1:1:T])
    highlight(h,[T+1:1:N],'NodeColor','r')
    
    filename = "data\"+num2str(T)+"_"+num2str(D)+"_"+num2str(K)+"_"+num2str(L)+"_"+"_Robot"+num2str(robot);
    saveas(h, filename,'jpeg')
end





