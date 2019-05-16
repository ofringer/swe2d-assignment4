%
% Set all computational edges to marker of type 0 and all other
% edges to 1, making sure to retain the type 2 or 3 edges
% previously specified in mark{u,v}face.
%
function boundary_edges %[markuface,markvface]=boundary_edges(cellmark,markuface,markvface)

global_pointers;

[Ni,Nj]=size(cellmark);
markuface0=markuface;
markvface0=markvface;

markuface=ones(Ni+1,Nj);
markvface=ones(Ni,Nj+1);
[is,js]=find(cellmark~=0); % Fixed on 5/28 - used to be
                           % find(cellmark~=1) but this missed the
                           % type 3 boundary on the right edge...
for m=1:length(is)
    i=is(m);
    j=js(m);
    if(i>1)
        if(cellmark(i-1,j)~=0)
            markuface(i,j)=0;
        end     
    end
    if(j>1)
        if(cellmark(i,j-1)~=0)
            markvface(i,j)=0;
        end
    end
end

markuface(markuface0==2)=2;
markuface(markuface0==3)=3;
markvface(markvface0==2)=2;
markvface(markvface0==3)=3;


