

#ifndef __ALTRACKHITEVENT_H__
#define __ALTRACKHITEVENT_H__

#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector2.h"




//------------------------------------------------------------------------------------
class AlTPCCluster : public TObject
{
    // TPC cluster data container
private:
    // Cluster properties
    Int_t track_id, row, sector;
    Float_t x_pos, y_pos, time;

public:
    AlTPCCluster() :
        track_id(0), row(0), sector(0), x_pos(0), y_pos(0), time(0)
    {

    }
        ~AlTPCCluster()
        {

        }

        // setters
        void set_cluster_pos_time(Float_t x, Float_t y, Float_t t) {x_pos = x; y_pos = y; time = t;}
        void set_track_id(Int_t i)                    {track_id = i; }
        void set_sector(Int_t i)                      {sector = i;   }
        void set_row(Int_t i)                         {row = i;      }

        // getters
        Float_t        get_cluster_x() const          { return x_pos;}
        Float_t        get_cluster_y() const          { return y_pos;}
        Float_t        get_cluster_time() const       { return time;}
        Int_t          get_cluster_row() const        { return row;}
        Int_t          get_cluster_sector() const     { return sector;}
        Int_t          get_cluster_track_id() const   { return track_id;}


        ClassDef(AlTPCCluster,1)  //
};
//------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------
class AlITSHit : public TObject
{
    // ITS hit data container
private:
    // ITS hit properties
    Float_t x_pos, y_pos, z_pos;
    int64_t BC;
    UShort_t row, col, id;

public:
    AlITSHit() :
        x_pos(0), y_pos(0), z_pos(0), BC(0), row(0), col(0), id(0)
    {

    }
        ~AlITSHit()
        {

        }

        // setters
        void set_cluster_pos(Float_t x, Float_t y, Float_t z) {x_pos = x; y_pos = y; z_pos = z;}
        void set_BC(int64_t i) {BC = i;}
        void set_row_col_id(UShort_t r, UShort_t c, UShort_t i) {row = r; col = c; id = i;}

        // getters
        Float_t        get_cluster_x() const          { return x_pos;}
        Float_t        get_cluster_y() const          { return y_pos;}
        Float_t        get_cluster_z() const          { return z_pos;}
        int64_t        get_BC()        const          { return BC;   }
        UShort_t       get_row()       const          { return row;  }
        UShort_t       get_col()       const          { return col;  }
        UShort_t       get_id()        const          { return id;   }


        ClassDef(AlITSHit,1)  //
};
//------------------------------------------------------------------------------------




//------------------------------------------------------------------------------------
class AlTrack : public TObject
{
private:
    // Track properties
    Int_t track_id;
    Float_t X, alpha, par[5], time0;
    Int_t fNumTPCCluster;

    TClonesArray* fTPCCluster; //->

public:
    AlTrack() :
        track_id(0), X(0), alpha(0), par(), time0(0), fNumTPCCluster(0)
    {
        fTPCCluster    = new TClonesArray( "AlTPCCluster", 10 );
    }
        ~AlTrack()
        {
            delete fTPCCluster;
            fTPCCluster = NULL;
        }

        // setters
        void set_track_id(Int_t i)               {track_id = i; }
        void set_X(Float_t f)                    {X = f; }
        void set_alpha(Float_t f)                {alpha = f; }
        void set_par(Float_t f0,Float_t f1,Float_t f2,Float_t f3,Float_t f4) {par[0] = f0; par[1] = f1; par[2] = f2; par[3] = f3; par[4] = f4;}
        void set_time0(Float_t f)                {time0 = f; }


        // getters
        Int_t    get_track_id() const        { return track_id;}
        Float_t  get_X() const               { return X;}
        Float_t  get_alpha() const           { return alpha;}
        Float_t  get_par(Int_t index) const  { return par[index];}
        Float_t  get_time0() const           { return time0;}


        //--------------------------------------
        // TPCCluster
        AlTPCCluster* createTPCCluster()
        {
            if (fNumTPCCluster == fTPCCluster->GetSize())
                fTPCCluster->Expand( fNumTPCCluster + 5 );
            if (fNumTPCCluster >= 5000000)
            {
                Fatal( "AlTPCCluster::createTPCCluster()", "ERROR: Too many TPCClusters (>5000000)!" );
                exit( 2 );
            }

            new((*fTPCCluster)[fNumTPCCluster++]) AlTPCCluster;
            return (AlTPCCluster*)((*fTPCCluster)[fNumTPCCluster - 1]);
        }
        void clearTPCClusterList()
        {
            fNumTPCCluster   = 0;
            fTPCCluster      ->Clear();
        }
        Int_t getNumTPCCluster() const
        {
            return fNumTPCCluster;
        }
        AlTPCCluster* getTPCCluster(Int_t i) const
        {
            return i < fNumTPCCluster ? (AlTPCCluster*)((*fTPCCluster)[i]) : NULL;
        }
        //--------------------------------------

        ClassDef(AlTrack,1)  // A simple track of a particle
};
//------------------------------------------------------------------------------------




//------------------------------------------------------------------------------------
class AlTrackHitEvent : public TObject
{
private:
    UShort_t      fNumTPCTrack;
    Int_t         fNumTPCCluster;
    Int_t         fNumITSHit;

    TClonesArray* fTPCTrack; //->
    TClonesArray* fTPCCluster; //->
    TClonesArray* fITSHit; //->

public:
    AlTrackHitEvent() :
        fNumTPCTrack(0),fNumTPCCluster(0),fNumITSHit(0)
    {
        fTPCTrack      = new TClonesArray( "AlTrack", 10 );
        fTPCCluster    = new TClonesArray( "AlTPCCluster", 10 );
        fITSHit        = new TClonesArray( "AlITSHit", 10 );
    }
        ~AlTrackHitEvent()
        {
            delete fTPCTrack;
            fTPCTrack = NULL;

            delete fTPCCluster;
            fTPCCluster = NULL;

            delete fITSHit;
            fITSHit = NULL;
        }

        // TPCTrack
        AlTrack* createTPCTrack()
        {
            if (fNumTPCTrack == fTPCTrack->GetSize())
                fTPCTrack->Expand( fNumTPCTrack + 5 );
            if (fNumTPCTrack >= 50000)
            {
                Fatal( "AlTrack::createTPCTrack()", "ERROR: Too many TPCTracks (>5000)!" );
                exit( 2 );
            }

            new((*fTPCTrack)[fNumTPCTrack++]) AlTrack;
            return (AlTrack*)((*fTPCTrack)[fNumTPCTrack - 1]);
        }
        void clearTPCTrackList()
        {
            fNumTPCTrack   = 0;
            fTPCTrack      ->Clear();
        }
        UShort_t getNumTPCTrack() const
        {
            return fNumTPCTrack;
        }
        AlTrack* getTPCTrack(UShort_t i) const
        {
            return i < fNumTPCTrack ? (AlTrack*)((*fTPCTrack)[i]) : NULL;
        }
        //--------------------------------------



        //--------------------------------------
        // TPCCluster
        AlTPCCluster* createTPCCluster()
        {
            if (fNumTPCCluster == fTPCCluster->GetSize())
                fTPCCluster->Expand( fNumTPCCluster + 5 );
            if (fNumTPCCluster >= 5000000)
            {
                Fatal( "AlTPCCluster::createTPCCluster()", "ERROR: Too many TPCClusters (>5000000)!" );
                exit( 2 );
            }

            new((*fTPCCluster)[fNumTPCCluster++]) AlTPCCluster;
            return (AlTPCCluster*)((*fTPCCluster)[fNumTPCCluster - 1]);
        }
        void clearTPCClusterList()
        {
            fNumTPCCluster   = 0;
            fTPCCluster      ->Clear();
        }
        Int_t getNumTPCCluster() const
        {
            return fNumTPCCluster;
        }
        AlTPCCluster* getTPCCluster(Int_t i) const
        {
            return i < fNumTPCCluster ? (AlTPCCluster*)((*fTPCCluster)[i]) : NULL;
        }
        //--------------------------------------


        //--------------------------------------
        // ITSHit
        AlITSHit* createITSHit()
        {
            if (fNumITSHit == fITSHit->GetSize())
                fITSHit->Expand( fNumITSHit + 5 );
            if (fNumITSHit >= 5000000)
            {
                Fatal( "AlITSHit::createITSHit()", "ERROR: Too many ITSHits (>5000000)!" );
                exit( 2 );
            }

            new((*fITSHit)[fNumITSHit++]) AlITSHit;
            return (AlITSHit*)((*fITSHit)[fNumITSHit - 1]);
        }
        void clearITSHitList()
        {
            fNumITSHit   = 0;
            fITSHit      ->Clear();
        }
        Int_t getNumITSHit() const
        {
            return fNumITSHit;
        }
        AlITSHit* getITSHit(Int_t i) const
        {
            return i < fNumITSHit ? (AlITSHit*)((*fITSHit)[i]) : NULL;
        }
        //--------------------------------------


        ClassDef(AlTrackHitEvent,1)  // A simple event compiled of tracks
};
//------------------------------------------------------------------------------------




#endif // __ALTRACKHITEVENT_H__