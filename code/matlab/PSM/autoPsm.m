function [psm,rp] = autoPsm(bimg,varargin)

% Automatically generates an ink-blob PSM for the input image.

% p = struct('skt',12,'spacing',4,'variance','proportional','sep',1);
% p = struct('skt',12,'spacing',4,'variance',0.5,'sep',1);

p = struct('skel',20,'spacing',4,'variance',0.5,'sep',1);

p = setupParameters(p,varargin{:});

% Check edge pixels to see if we need to invert
epix = [bimg(1,:),bimg(end,:),bimg(2:end-1,1)',bimg(2:end-1,end)'];
if (nnz(epix)>numel(epix)/2)&&~isfield(p,'sic')
    warning('Inverting input image');
    bimg = ~bimg;
end;

if isempty(bimg)
    warning ('PSM of empty image');
    bimg = zeros(100);
end;

[nrow,ncol] = size(bimg);
if isfield(p,'skel')&&p.skel
    sk = bimg;
    skr = sk;  % radius is unknown
else
    %[skg,skr] = skelgrad(bimg);  % has been crashing...
    skr = double(bwdist(imdilate(~bimg,strel(eye(3))))).^2;  % exact
    skg = 2*p.skt*bwmorph(bimg,'thin',inf);  % approximation

    
    skr = bimg+sqrt(skr);  % this gives the actual distance to the edge
    sk = bwmorph(skg>p.skt,'thin',inf);
    if all(sk(:)==0)
        sk = bwmorph(skg==max(skg(:)),'thin',inf);
    end;
end;

% analyze
psk = padarray(sk,[1 1],0);
nnbr = psk(1:end-2,1:end-2)+psk(1:end-2,2:end-1)+psk(1:end-2,3:end)...
    +psk(2:end-1,1:end-2)+psk(2:end-1,3:end)...
    +psk(3:end,1:end-2)+psk(3:end,2:end-1)+psk(3:end,3:end);
sknbr = nnbr(sk);
eid = find(sk&nnbr<2)';   % ends
bid = find(sk&nnbr==2)';  % branch lines
%jid = find(sk&nnbr>2)';   % junctions -- replaced with below

% work on junctions to pick just one from each connected component
[jcc,njcc] = bwlabel(sk&nnbr>2);
jid = zeros(1,njcc);
for j = 1:njcc
    jjid = find(jcc==j);
    [ji,jj] = ind2sub(size(sk),jjid);
    jd = (ji-mean(ji)).^2+(jj-mean(jj)).^2;
    [~,jdmi] = min(jd);
    jid(j) = jjid(jdmi);
end;

% uniform spacing replaces radius everywhere with median value
if isfield(p,'spacing')&&~isempty(p.spacing)
    skr0 = skr;
    if ischar(p.spacing)&&strcmp(p.spacing,'uniform')
        %skr(:) = median(skr(sk));
        skr(:) = trimmean(skr(sk),50);
    else
        skr(:) = p.spacing;
    end;
end;

% block out end/junction points
pick = [eid jid];
d = inf(nrow,ncol);
src = nan(nrow,ncol);
for i = 1:numel(pick)
    j = pick(i);
    m = repmat(1e10,nrow,ncol);
    m(j) = 0;
    m = gdt2d(m,skr(j).^2,skr(j).^2);
    src(m<=d) = j;
    d = min(d,m);
end;

% add more points
skd = sk.*d;
th = 2*p.sep;
while any(skd(:)>=th)
    j = find(skd==min(skd(skd>=th)),1);
    pick(end+1) = j;
    m = repmat(1e10,nrow,ncol);
    m(j) = 0;
    m = gdt2d(m,skr(j).^2,skr(j).^2);
    src(m<=d) = j;
    d = min(d,m);
    skd = sk.*d;
    if ~any(skd(:)>=th)&&(~isfield(p,'stragglers')||p.stragglers)
        % lower threshold at end to get stragglers
        th = p.sep;
    end;
end;

% eliminate points outside mask, if provided
if isfield(p,'mask')
    pick = pick(mask(pick));
end;

% now we have all our picks -- build tree based on distance
npick = numel(pick);
pickp = zeros(nrow,ncol);
pickp(pick) = 1;
[ip,jp] = find(pickp);
psm = struct('x',0,'y',0,'r',num2cell(skr(sort(pick))),...
    'parent',0,'children',[],'absx',0,'absy',0,...
    'vx',1,'vy',1);
if ~isfield(p,'root')
    cip = mean(ip);
    cjp = mean(jp);
else
    cip = p.root(2);
    cjp = p.root(1);
end;
sqcdp = (ip-cip).^2+(jp-cjp).^2;
[~,root] = min(sqcdp);
sqpd = (ip-ip(root)).^2+(jp-jp(root)).^2;
src = repmat(root,1,npick);
treed = false(npick,1);
psm(root) = struct('x',jp(root)-cjp,'y',ip(root)-cip,...
    'r',psm(root).r,'parent',0,'children',[],...
    'absx',jp(root)-cjp,'absy',ip(root)-cip,...
    'vx',1,'vy',1);
treed(root) = true;
while ~all(treed)
    [~,node] = min(sqpd(:)+1e10*treed);
    parent = src(node);
    sqipd = (ip-ip(node)).^2+(jp-jp(node)).^2;
    src(sqipd(:)<sqpd(:)) = node;
    sqpd = min(sqpd(:),sqipd(:));
    treed(node) = true;
    psm(node) = struct('x',jp(node)-jp(parent),'y',ip(node)-ip(parent),...
        'r',psm(node).r,'parent',parent,'children',[],...
        'absx',jp(node)-cjp,'absy',ip(node)-cip,...
        'vx',1,'vy',1);
    psm(parent).children(end+1) = node;
end;
rp = [jp(root);ip(root)];

if isfield(p,'radrestore')&&p.radrestore
    radii = num2cell(skr0(ip+nrow.*jp-nrow));
    [psm.r] = deal(radii{:});
end;
if isfield(p,'radius')
    [psm.r] = deal(p.radius);
end;

% proportional variance is relative to parent distance
if isfield(p,'variance')&&ischar(p.variance)&&strcmp(p.spacing,'proportional')
    x = [psm.x];
    y = [psm.y];
    d = sqrt(x.^2+y.^2);
    pv = num2cell(d./mean(d));
    [psm.vx] = deal(pv{:});
    [psm.vy] = deal(pv{:});
end;

% replace variance with specified value
if isfield(p,'variance')&&isnumeric(p.variance)
    [psm.vx] = deal(p.variance);
    [psm.vy] = deal(p.variance);
end;

% experimental:  weaken ties between separate connected components
if isfield(p,'ccweaken')
    lbl = bwlabel(bimg);
    plbl = lbl(pick);
    for i = 1:numel(psm)
        parent = psm(i).parent;
        if (parent > 0)&&(plbl(i)~=plbl(parent))
            psm(i).vx = psm(i).vx*p.ccweaken;
            psm(i).vy = psm(i).vy*p.ccweaken;
        end;
    end;
end;

if isfield(p,'code')&&(islogical(p.code)||isnumeric(p.code))&&p.code
    spick = sort(pick);
    [psm.code] = deal(1);
    [psm(ismember(spick,eid)).code] = deal(2);
    [psm(ismember(spick,jid)).code] = deal(3);
elseif isfield(p,'code')&&ischar(p.code)&&strcmp(p.code,'descript')
    sd = skelDescriptor(sk);
    [psm.code] = celldeal(num2cell(sd(sort(pick))));
elseif isfield(p,'code')&&ischar(p.code)&&strcmp(p.code,'theta')
    sd = skelDescriptor(sk);
    [psm.code] = celldeal(num2cell(sd(sort(pick))));
end;

if isfield(p,'rsig')
    [psm.rsig] = deal(p.rsig);
end;
end


function [varargout] = celldeal(c)
[varargout{1:numel(c)}] = deal(c{:});
end


function E = even(M)
E = (round(M/2) == (M/2));
end


function y = linsmooth(x,sigma,endmethod)
if (nargin < 2)
    sigma = 1;
end;
if (nargin < 3)
    endmethod = 'replicate';
end;

if (sigma == 0)
    y = x;
else
    radius = ceil(2.5*sigma);
    h = fspecial('gaussian',[2*radius+1,1],sigma);
    px = padarray(shiftdim(x),radius,endmethod,'both');
    yc = conv(px,h);
    y = reshape(yc(2*radius+1:end-2*radius),size(x));
end;
end


function p = setupParameters(defp,varargin)

if (nargout == 0)
    warning('setupParameters called with no outputs.');
end;
p = defp;
if even(nargin)
    p2 = varargin{1};
    fn = fieldnames(p2);
    for i = 1:numel(fn)
        p.(fn{i}) = p2.(fn{i});
    end;
elseif isfield(defp,'default')
    p = 'default';
end;
if ischar(p)
    p = defp.(p);
end;
for i = (1+even(nargin)):2:nargin-2
    p(1).(varargin{i}) = varargin{i+1};
end;
end


% Returns coordinates of each 1D section of a skeleton image.  
% Result is a cell vector of 2xP xy coordinates.
function [arcs,joints,ends] = skelArcs(tbimg)
% analyze
[nrow,ncol] = size(tbimg);
pimg = padarray(tbimg,[1 1],0);
nnbr = pimg(1:end-2,1:end-2)+pimg(1:end-2,2:end-1)+pimg(1:end-2,3:end)...
    +pimg(2:end-1,1:end-2)+pimg(2:end-1,3:end)...
    +pimg(3:end,1:end-2)+pimg(3:end,2:end-1)+pimg(3:end,3:end);
sknbr = nnbr(logical(tbimg));

pnbr = padarray(nnbr,[1 1],0);
peid = find(pimg&(pnbr<2))';   % ends
pbid = find(pimg&(pnbr==2))';  % branch lines
pjid = find(pimg&(pnbr>2))';   % junctions -- may be small connected component

if (numel(pjid)>numel(pbid))
    warning('Too many junctions -- input image not thinned?');
end;

arcs = cell(1,sum(pnbr([pjid peid])));
narc = 0;
unused = pimg&padarray(nnbr<=2,[1 1],false);
[xg1,yg1] = meshgrid(0:ncol+1,0:nrow+1);
[xg,yg] = meshgrid(1:ncol,1:nrow);
ig1 = reshape(1:(nrow+2)*(ncol+2),nrow+2,ncol+2);
nhood = bsxfun(@plus,[-nrow-2,0,nrow+2],[-1 0 1]');
jmap = zeros(nrow+2,ncol+2);
jmap(pjid) = 1;
for i = [pjid peid pbid]
    jmap(i) = 0;
    if unused(i)  % necessary if we start on a loop
        unused(i) = false;
    end;
    local = i+nhood;
    while any(unused(local(:)))
        j = ig1(i+nhood(find(unused(local(:)),1)));
        unused(j) = false;
        narc = narc+1;
        arcs{narc} = zeros(1,nrow*ncol);
        arcs{narc}(1) = i;
        arcs{narc}(2) = j;
        arclen = 2;
        local = j+nhood;
        while any(unused(local(:)))
            j = ig1(j+nhood(find(unused(local(:)),1)));
            unused(j) = false;
            arclen = arclen+1;
            arcs{narc}(arclen) = j;
            local = j+nhood;
        end;
        if any(jmap(local(:))) % add final junction, if any
            j = ig1(j+nhood(find(jmap(local(:)),1)));
            arclen = arclen+1;
            arcs{narc}(arclen) = j;
        end;
        arcs{narc} = [xg1(arcs{narc}(1:arclen));yg1(arcs{narc}(1:arclen))];
        local = i+nhood;
    end;
    jmap(pjid) = 1;
end;
arcs = arcs(1:narc);

if (nargout > 1)
    % produce list of joints as well
    jrp = regionprops(tbimg&(nnbr>2),'PixelList');
    joints  = cellfun(@(a)a',{jrp.PixelList},'uniform',0);
end;

if (nargout > 2)
    % produce list of ends as well
    [erpi,erpj] = find(tbimg&(nnbr<=1));
    ends = num2cell([erpj';erpi'],1);
    %erp = regionprops(nnbr<=1,'PixelList');
    %ends  = cellfun(@(a)a',{erp.PixelList},'uniform',0);
end;
end


% Return descriptor codes for the supplied skeleton image.
%   0 = Not on skeleton
%   1 = Junction
%   2 = Endpoint
%   pi..2*pi = Arc angle
function sd = skelDescriptor(tbimg,thsig,markJoints,markEnds)

if (nargin < 2)||isempty(thsig)
    thsig = 2;
end;
if (nargin < 3)||isempty(markJoints)
    markJoints = false;
end;
if (nargin < 4)||isempty(markEnds)
    markEnds = false;
end;

% Get arc, joint, end point lists
[arcs,joints,ends] = skelArcs(tbimg);

sd = zeros(size(tbimg));

% Compute & smooth arc curvature
narc = numel(arcs);
ath = cell(1,narc);
for i = 1:narc
    ath{i} = atan2(diff(arcs{i}(2,:)),diff(arcs{i}(1,:)));
    athc = zeros(1,numel(ath{i}));
    athc(find(diff(ath{i})>pi)+1) = -2*pi;
    athc(find(diff(ath{i})<-pi)+1) = 2*pi;
    ath{i} = ath{i}+cumsum(athc);
    ath{i} = conv(ath{i}([1 1:end end]),[0.5 0.5],'valid');
    ath{i} = linsmooth(ath{i},thsig);
    sd(arcs{i}(2,:)+size(sd,1)*(arcs{i}(1,:)-1)) = mod(ath{i},pi)+pi;
end;

if markJoints
    j = [joints{:}];
    sd(j(2,:)+size(sd,1)*(j(1,:)-1)) = 1;
end;

if markEnds
    e = [ends{:}];
    sd(e(2,:)+size(sd,1)*(e(1,:)-1)) = 2;
end;
end
