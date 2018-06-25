%% Visualization of the path generated for each robot

function [] = map(A, robot, T, N, x_pos, y_pos)
    figure;
    G = digraph(A);
    h = plot(G,'XData', x_pos, 'YData', y_pos);
    axis square
    hold on
    title("Path of Robot "+num2str(robot));
    highlight(h,[1:1:T])
    highlight(h,[T+1:1:N],'NodeColor','r')
end





