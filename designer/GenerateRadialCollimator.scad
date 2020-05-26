include<RadialCollimatorGrid.scad>;  //Defines the variable mywalls and
include<RadialCollimatorBox.scad>;  //Defines the variable mybox


halfWallThickness=.125;  // in mm

BoxWallTh=.5;  // in mm
    
function faceNormal(v1,v2) = cross(v1,v2)/norm(cross(v1,v2));

module drawWall(l)
{
    fn=faceNormal(l[0]-l[1],l[0]-l[2]);
    fsw=halfWallThickness*fn;
    points= [l[3]-fsw,l[0]-fsw,l[1]-fsw,l[2]-fsw,l[3]+fsw,l[0]+fsw,l[1]+fsw,l[2]+fsw];
    faces=[[0,1,2,3],[4,5,1,0],[7,6,5,4],[5,6,2,1],[6,7,3,2],[7,4,0,3]];
    polyhedron ( points,faces );
}

module drawStructure(){
for(a=mywalls)
    {
        drawWall(a);
    }

}


// Clipping...
pointsClip= [c2ThetaMinPhiMax,c1ThetaMinPhiMax,c1ThetaMaxPhiMax,c2ThetaMaxPhiMax,c2ThetaMinPhiMin,c1ThetaMinPhiMin,c1ThetaMaxPhiMin,c2ThetaMaxPhiMin];
facesClip=[[0,1,2,3],[4,5,1,0],[7,6,5,4],[5,6,2,1],[6,7,3,2],[7,4,0,3]];
intersection(){
drawStructure();
polyhedron ( pointsClip,facesClip );
}


module drawBoxWall(p)
{
    fn=p[1];
    l=p[0];
    fsw=BoxWallTh*fn;
    points= [l[3],l[0],l[1],l[2],l[3]+fsw,l[0]+fsw,l[1]+fsw,l[2]+fsw];
    faces=[[0,1,2,3],[4,5,1,0],[7,6,5,4],[5,6,2,1],[6,7,3,2],[7,4,0,3]];
    polyhedron ( points,faces );
}

drawBoxWall(mybox[0]);
drawBoxWall(mybox[1]);
drawBoxWall(mybox[2]);
drawBoxWall(mybox[3]);


echo(c2ThetaMinPhiMax);
