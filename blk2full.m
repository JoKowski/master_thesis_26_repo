function M = blk2full(B)
%BLK2FULL  N x N Blockmatrix (Cell mit 2x2-Bloecken) -> ausgeschriebene Matrix.
%
%   M = BLK2FULL(B)   B: N x N cell, jede Zelle 2x2 (oder leer = Nullblock)
%                     M: 2N x 2N
%
%   Y_b ist nach Gu et al. Gl.(15) eine N x N BLOCKmatrix; jeder Eintrag
%   Y_bkl ist die 2x2-Matrix aus Gl.(16).
%
%   Ausgeschrieben wird nur dort, wo gerechnet werden muss: inv, det, eig,
%   lu und mldivide arbeiten auf der 2N x 2N-Matrix. Blockweise ginge das
%   nur, wenn die Bloecke untereinander kommutieren - sie tun es nicht.
%
%   Siehe auch FULL2BLK.

validateattributes(B, {'cell'}, {'2d','square'}, mfilename, 'B', 1);
N = size(B,1);
proto = [];
for k = 1:N
    for l = 1:N
        if ~isempty(B{k,l}), proto = B{k,l}; break; end
    end
    if ~isempty(proto), break; end
end
if isempty(proto), proto = zeros(2,2); end
if isa(proto,'sym')
    M = sym(zeros(2*N,2*N));
else
    M = zeros(2*N,2*N); 
end

for k = 1:N
    for l = 1:N
        if isempty(B{k,l}), continue; end
        validateattributes(B{k,l}, {'numeric','sym'}, {'size',[2 2]}, ...
                           mfilename, sprintf('B{%d,%d}',k,l));
        M(2*k-1:2*k, 2*l-1:2*l) = B{k,l};
    end
end
end
