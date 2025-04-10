close all;
%%
figure
for k = 1:400
    plot(p1000.Bg_coupleSeq{1,k,1}), hold on
end

%%
figure
for k = 1:400
    plot(p1000.Bg_selfSeq{1,k,1}), hold on
end

%%
figure
for k = 1:100
    plot(Bg_coupleSeq{1,k,3}), hold on
end

%%
figure
for k = 1:12
    plot(Bg(5).field{1,k}), hold on
end


%%
figure
subplot(3,1,1)
for k = 1:50
    plot(Bg_Seq{1,k,1}), hold on
end
subplot(3,1,2)
for k = 1:50
    plot(BgSeq{1,k,2}), hold on
end
subplot(3,1,3)
for k = 1:50
    plot(BgSeq{1,k,3}), hold on
end

%%
figure
for k = 1:50
    plot(Bg{2,k}-Bg2{2,k}), hold on
end

%%
figure
for k = 1:50
    plot(Bg_couple{2,k}-Bg2{2,k}), hold on
end

%%
figure
for k = 1:50
    plot(Bg_selfSeq{1,k,1}-BG_self(k,1).shape), hold on
end

%%
figure
for k = 1:50
    plot(Bg_selfSeq{2,k,1}-BG_coup(k,1).shape), hold on
end

%%
figure
for k = 1:50
    plot(Bg_coupleSeq{1,k,1}-BG_self(k,1).shape), hold on
end

%%
figure
subplot(3,1,1)
for k = 1:50
    plot(Bg_coupleSeq{2,k,1}), hold on
end
subplot(3,1,2)
for k = 1:50
    plot(Bg_coupleSeq{2,k,2}), hold on
end
subplot(3,1,3)
for k = 1:50
    plot(Bg_coupleSeq{2,k,3}), hold on
end

%%
figure
for k = 1:12
    plot(Bg(5).field{k}), hold on
end

%%
figure
for k = 1:1:400
    plot(p1000.Bg_selfSeq{2,k,3}), hold on
end

%% 
phase_sum = 0;
for k =1:400
    phase_sum =phase_sum + exp(1i*p1000.Bg_selfSeq{1,k,3}(500));
end
disp(phase_sum);
