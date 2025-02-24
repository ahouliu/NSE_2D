%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Poincare map Ek(tn) when dEk(tn)=0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EkP,err] = Pmap_func(Ek,dEk)

    if length(Ek)~=length(dEk)
        error('not match');
    end

    LEN  = length(Ek);
    EkP  = zeros(LEN,1);
    dEkP = zeros(LEN,1);

    for ind=1:LEN-1
        % up zero-cross
        if ( dEk(ind)<=0 && dEk(ind+1)>0 ) || ( dEk(ind)<0 && dEk(ind+1)>=0 )
            loc1 = ind;
            loc2 = ind+1;
            loc0 = loc1-dEk(loc1)/(dEk(loc2)-dEk(loc1));
            EkP(ind) = Ek(loc1) - (loc1-loc0)*(Ek(loc2)-Ek(loc1));
            continue;
        end
        % down zero-corss
        if ( dEk(ind)>=0 && dEk(ind+1)<0 ) || ( dEk(ind)>0 && dEk(ind+1)<=0 )
            loc1 = ind;
            loc2 = ind+1;
            loc0 = loc1-dEk(loc1)/(dEk(loc2)-dEk(loc1));
            EkP(ind) = Ek(loc1) - (loc1-loc0)*(Ek(loc2)-Ek(loc1));
            continue;
        end  

        % zero
        if abs(dEk(ind))<1e-10
            EkP(ind)  = ind;
            dEkP(ind) = dEk(ind);
            continue;
        end
        EkP(ind) = NaN;

    end
    EkP(isnan(EkP)) = [];
    err = max(abs(dEkP));

end