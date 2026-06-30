function lz = lzc_norm(x, binarize)

if nargin < 2 || isempty(binarize), binarize = 'median'; end

x = double(x(:));
x = x(isfinite(x));
n = numel(x);

if n < 20
    lz = NaN;
    return;
end

switch lower(binarize)
    case 'median'
        thr = median(x);
    case 'mean'
        thr = mean(x);
    case 'zero'
        thr = 0;
    otherwise
        thr = median(x);
end

s = x > thr;
s = uint8(s(:)');   % row of 0/1

% --- LZ76 complexity ---
c = 1;
i = 1; k = 1; l = 1;

while true
    if i + k > n
        c = c + 1;
        break;
    end

    sub1 = s(i:i+k-1);
    sub2 = s(l:l+k-1);

    if isequal(sub1, sub2)
        k = k + 1;
        if l + k - 1 > n
            c = c + 1;
            break;
        end
    else
        l = l + 1;
        if l == i
            c = c + 1;
            i = i + k;
            if i > n
                break;
            end
            l = 1;
            k = 1;
        end
    end
end

% --- normalization (common choice) ---
% b(n) = n/log2(n)  (for binary alphabet)
bn = n / log2(n);
lz = c / bn;
end
