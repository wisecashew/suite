#!/bin/awk -f

BEGIN {prior=1; atoms_section=0; velocities_section=0; bonds_section=0; angles_section=0; post=0}

/^$/                       {print; next;}
/^Atoms/                   {prior=0; atoms_section=1; print; next;}
/^Velocities/              {atoms_section=0; velocities_section=1; print; next;}
/^Bonds/                   {velocities_section=0; bonds_section=1; print; next;}
/^Angles/                  {bonds_section=0; angles_section=1; print; next;}
/^Dihedrals/               {angles_section=0; post=1; print; next}
prior                      {print; next}
atoms_section  && $2 == 1  {print; next;}
bonds_section  && $2 < 15  {print; next;}
angles_section && $2 < 31  {print; next;}
post                       {print; next;}

END {}
