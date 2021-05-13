-- module_adam_feedbf_les_press_velfg_ve_etc_superkernel
module ASTInstance ( ast
        , functionSignaturesList
        , stencilDefinitionsList
        , mainArgDeclsList 
        , scalarisedArgsList
        , origNamesList
        ) where
import TyTraCLAST

ast :: TyTraCLAST
ast = [
        -- adam_map_36
       ( (Tuple [Vec VT (Scalar VDC DFloat "f_1"),Vec VT (Scalar VDC DFloat "fold_1"),Vec VT (Scalar VDC DFloat "g_1"),Vec VT (Scalar VDC DFloat "gold_1"),Vec VT (Scalar VDC DFloat "h_1"),Vec VT (Scalar VDC DFloat "hold_1")]), UnzipT ( Map (Function "adam_map_36" []) (ZipT [Vec VI (Scalar VDC DFloat "f_0"),Vec VI (Scalar VDC DFloat "g_0"),Vec VI (Scalar VDC DFloat "h_0"),Vec VI (Scalar VDC DFloat "fold_0"),Vec VI (Scalar VDC DFloat "gold_0"),Vec VI (Scalar VDC DFloat "hold_0")]) ) )
        -- feedbf_map_49
       ,( (Tuple [Vec VO (Scalar VDC DFloat "fx_0"),Vec VO (Scalar VDC DFloat "fy_0"),Vec VO (Scalar VDC DFloat "fz_0"),Vec VT (Scalar VDC DFloat "usum_1"),Vec VT (Scalar VDC DFloat "vsum_1"),Vec VT (Scalar VDC DFloat "wsum_1")]), UnzipT ( Map (Function "feedbf_map_49"  [Scalar VI DFloat "alpha_0",Scalar VI DFloat "dt_0",Scalar VI DFloat "beta_0"]) (ZipT [Vec VI (Scalar VDC DFloat "usum_0"),Vec VI (Scalar VDC DFloat "u_0"),Vec VI (Scalar VDC DFloat "bmask1_0"),Vec VI (Scalar VDC DFloat "vsum_0"),Vec VI (Scalar VDC DFloat "v_0"),Vec VI (Scalar VDC DFloat "cmask1_0"),Vec VI (Scalar VDC DFloat "wsum_0"),Vec VI (Scalar VDC DFloat "w_0"),Vec VI (Scalar VDC DFloat "dmask1_0")]) ) )
        -- feedbf_map_67
       ,( (Tuple [Vec VT (Scalar VDC DFloat "f_2"),Vec VT (Scalar VDC DFloat "g_2"),Vec VT (Scalar VDC DFloat "h_2")]), UnzipT ( Map (Function "feedbf_map_67" []) (ZipT [Vec VT (Scalar VDC DFloat "f_1"),Vec VO (Scalar VDC DFloat "fx_0"),Vec VT (Scalar VDC DFloat "g_1"),Vec VO (Scalar VDC DFloat "fy_0"),Vec VT (Scalar VDC DFloat "h_1"),Vec VO (Scalar VDC DFloat "fz_0")]) ) )
        -- les_map_119
       ,( Vec VO (Scalar VDC DFloat "sm_0"), Map (Function "les_map_119"  [Scalar VI DFloat "dx1_0",Scalar VI DFloat "dy1_0",Scalar VI DFloat "dzs_0",Scalar VI DFloat "dzn_0",Scalar VI DFloat "delx1_0"]) (ZipT [Vec VI (Scalar VDC DFloat "u_0"),Vec VI (Scalar VDC DFloat "v_0"),Vec VI (Scalar VDC DFloat "w_0")]) )
        -- les_map_164
       ,( Vec VT (Scalar VDC DFloat "f_3"), Map (Function "les_map_164"  [Scalar VI DFloat "dy1_0",Scalar VI DFloat "dx1_0",Scalar VI DFloat "dzn_0",Scalar VI DFloat "dzs_0",Scalar VI DFloat "dxs_0"]) (ZipT [Vec VO (Scalar VDC DFloat "sm_0"),Vec VI (Scalar VDC DFloat "u_0"),Vec VI (Scalar VDC DFloat "v_0"),Vec VI (Scalar VDC DFloat "w_0"),Vec VT (Scalar VDC DFloat "f_2")]) )
        -- les_map_202
       ,( Vec VT (Scalar VDC DFloat "g_3"), Map (Function "les_map_202"  [Scalar VI DFloat "dy1_0",Scalar VI DFloat "dx1_0",Scalar VI DFloat "dzn_0",Scalar VI DFloat "dzs_0",Scalar VI DFloat "dys_0"]) (ZipT [Vec VO (Scalar VDC DFloat "sm_0"),Vec VI (Scalar VDC DFloat "u_0"),Vec VI (Scalar VDC DFloat "v_0"),Vec VI (Scalar VDC DFloat "w_0"),Vec VT (Scalar VDC DFloat "g_2")]) )
        -- les_map_240
       ,( Vec VT (Scalar VDC DFloat "h_3"), Map (Function "les_map_240"  [Scalar VI DFloat "dzn_0",Scalar VI DFloat "dx1_0",Scalar VI DFloat "dy1_0",Scalar VI DFloat "dzs_0"]) (ZipT [Vec VO (Scalar VDC DFloat "sm_0"),Vec VI (Scalar VDC DFloat "u_0"),Vec VI (Scalar VDC DFloat "w_0"),Vec VI (Scalar VDC DFloat "v_0"),Vec VT (Scalar VDC DFloat "h_2")]) )
        -- press_map_65
       ,( Vec VT (Scalar VDC DFloat "rhs_1"), Map (Function "press_map_65"  [Scalar VI DFloat "dx1_0",Scalar VI DFloat "dy1_0",Scalar VI DFloat "dzn_0",Scalar VI DFloat "dt_0"]) (ZipT [Vec VI (Scalar VDC DFloat "u_0"),Vec VI (Scalar VDC DFloat "v_0"),Vec VI (Scalar VDC DFloat "w_0"),Vec VT (Scalar VDC DFloat "f_3"),Vec VT (Scalar VDC DFloat "g_3"),Vec VT (Scalar VDC DFloat "h_3"),Vec VT (Scalar VDC DFloat "rhs_0")]) )
        -- press_reduce_78
       , ((Tuple [Scalar VO DFloat "area_1",Scalar VO DFloat "rhsav_1"]), Fold (Function "press_reduce_78"  [Scalar VI DFloat "dx1_0",Scalar VI DFloat "dy1_0",Scalar VI DFloat "dzn_0"]) (Scalar VI DFloat "rhsav_0",Scalar VI DFloat "area_0") (Vec VT (Scalar VDC DFloat "rhs_1")))
        -- press_map_87
       ,( Vec VO (Scalar VDC DFloat "rhs_2"), Map (Function "press_map_87"  [Scalar VI DFloat "rhsav_1"]) (Vec VT (Scalar VDC DFloat "rhs_1")) )
        -- press_map_97
       ,( (Tuple [Vec VT (Scalar VDC DFloat "p0_1"),Vec VO (Scalar VDC DFloat "p1_1")]), UnzipT ( Map (Function "press_map_97"  [Scalar VI DFloat "dzs_0",Scalar VI DFloat "dys_0",Scalar VI DFloat "dxs_0",Scalar VI DInt "nrd_0"]) (ZipT [Vec VI (Scalar VDC DFloat "p0_0"),Vec VO (Scalar VDC DFloat "rhs_2"),Vec VT (Scalar VDC DFloat "p1_0")]) ) )
        -- press_reduce_129
       , ((Tuple [Scalar VO DFloat "pav_1",Scalar VO DFloat "pco_1"]), Fold (Function "press_reduce_129"  [Scalar VI DFloat "dx1_0",Scalar VI DFloat "dy1_0",Scalar VI DFloat "dzn_0"]) (Scalar VI DFloat "pav_0",Scalar VI DFloat "pco_0") (Vec VT (Scalar VDC DFloat "p0_1")))
        -- press_map_138
       ,( Vec VT (Scalar VDC DFloat "p0_2"), Map (Function "press_map_138"  [Scalar VI DFloat "pav_1"]) (Vec VT (Scalar VDC DFloat "p0_1")) )
        -- velfg_map_95
       ,( Vec VO (Scalar VDC DFloat "f_4"), Map (Function "velfg_map_95"  [Scalar VI DFloat "dx1_0",Scalar VI DFloat "dy1_0",Scalar VI DFloat "dzs_0",Scalar VI DFloat "dzn_0"]) (ZipT [Vec VI (Scalar VDC DFloat "u_0"),Vec VI (Scalar VDC DFloat "v_0"),Vec VI (Scalar VDC DFloat "w_0")]) )
        -- velfg_map_135
       ,( Vec VO (Scalar VDC DFloat "g_4"), Map (Function "velfg_map_135"  [Scalar VI DFloat "dy1_0",Scalar VI DFloat "dx1_0",Scalar VI DFloat "dzs_0",Scalar VI DFloat "dzn_0"]) (ZipT [Vec VI (Scalar VDC DFloat "u_0"),Vec VI (Scalar VDC DFloat "v_0"),Vec VI (Scalar VDC DFloat "w_0")]) )
        -- velfg_map_175
       ,( Vec VO (Scalar VDC DFloat "h_4"), Map (Function "velfg_map_175"  [Scalar VI DFloat "dzn_0",Scalar VI DFloat "dx1_0",Scalar VI DFloat "dy1_0"]) (ZipT [Vec VI (Scalar VDC DFloat "u_0"),Vec VI (Scalar VDC DFloat "w_0"),Vec VI (Scalar VDC DFloat "v_0")]) )
        -- velnw_map_36
       ,( Vec VT (Scalar VDC DFloat "u_1"), Map (Function "velnw_map_36"  [Scalar VI DFloat "ro_0",Scalar VI DFloat "dxs_0",Scalar VI DFloat "dt_0"]) (ZipT [Vec VT (Scalar VDC DFloat "p0_2"),Vec VI (Scalar VDC DFloat "u_0"),Vec VO (Scalar VDC DFloat "f_4")]) )
        -- velnw_map_44
       ,( Vec VT (Scalar VDC DFloat "v_1"), Map (Function "velnw_map_44"  [Scalar VI DFloat "ro_0",Scalar VI DFloat "dys_0",Scalar VI DFloat "dt_0"]) (ZipT [Vec VT (Scalar VDC DFloat "p0_2"),Vec VI (Scalar VDC DFloat "v_0"),Vec VO (Scalar VDC DFloat "g_4")]) )
        -- velnw_map_52
       ,( Vec VT (Scalar VDC DFloat "w_1"), Map (Function "velnw_map_52"  [Scalar VI DFloat "ro_0",Scalar VI DFloat "dzs_0",Scalar VI DFloat "dt_0"]) (ZipT [Vec VT (Scalar VDC DFloat "p0_2"),Vec VI (Scalar VDC DFloat "w_0"),Vec VO (Scalar VDC DFloat "h_4")]) )
        ]

functionSignaturesList = [
        ("adam_map_36",  [Tuple [],Tuple [Scalar VDC DFloat "f_0",Scalar VDC DFloat "g_0",Scalar VDC DFloat "h_0",Scalar VDC DFloat "fold_0",Scalar VDC DFloat "gold_0",Scalar VDC DFloat "hold_0"],Tuple [Scalar VDC DFloat "f_1",Scalar VDC DFloat "fold_1",Scalar VDC DFloat "g_1",Scalar VDC DFloat "gold_1",Scalar VDC DFloat "h_1",Scalar VDC DFloat "hold_1"]]),
        ("feedbf_map_49",  [Tuple [Scalar VDC DFloat "alpha_0",Scalar VDC DFloat "dt_0",Scalar VDC DFloat "beta_0"],Tuple [Scalar VDC DFloat "usum_0",Scalar VDC DFloat "u_0",Scalar VDC DFloat "bmask1_0",Scalar VDC DFloat "vsum_0",Scalar VDC DFloat "v_0",Scalar VDC DFloat "cmask1_0",Scalar VDC DFloat "wsum_0",Scalar VDC DFloat "w_0",Scalar VDC DFloat "dmask1_0"],Tuple [Scalar VDC DFloat "fx_0",Scalar VDC DFloat "fy_0",Scalar VDC DFloat "fz_0",Scalar VDC DFloat "usum_1",Scalar VDC DFloat "vsum_1",Scalar VDC DFloat "wsum_1"]]),
        ("feedbf_map_67",  [Tuple [],Tuple [Scalar VDC DFloat "f_1",Scalar VDC DFloat "fx_0",Scalar VDC DFloat "g_1",Scalar VDC DFloat "fy_0",Scalar VDC DFloat "h_1",Scalar VDC DFloat "fz_0"],Tuple [Scalar VDC DFloat "f_2",Scalar VDC DFloat "g_2",Scalar VDC DFloat "h_2"]]),
        ("les_map_119",  [Tuple [Scalar VDC DFloat "dx1_0",Scalar VDC DFloat "dy1_0",Scalar VDC DFloat "dzs_0",Scalar VDC DFloat "dzn_0",Scalar VDC DFloat "delx1_0"],Tuple [Scalar VDC DFloat "u_0",Scalar VDC DFloat "v_0",Scalar VDC DFloat "w_0"],Scalar VDC DFloat "sm_0"]),
        ("les_map_164",  [Tuple [Scalar VDC DFloat "dy1_0",Scalar VDC DFloat "dx1_0",Scalar VDC DFloat "dzn_0",Scalar VDC DFloat "dzs_0",Scalar VDC DFloat "dxs_0"],Tuple [Scalar VDC DFloat "sm_0",Scalar VDC DFloat "u_0",Scalar VDC DFloat "v_0",Scalar VDC DFloat "w_0",Scalar VDC DFloat "f_2"],Scalar VDC DFloat "f_3"]),
        ("les_map_202",  [Tuple [Scalar VDC DFloat "dy1_0",Scalar VDC DFloat "dx1_0",Scalar VDC DFloat "dzn_0",Scalar VDC DFloat "dzs_0",Scalar VDC DFloat "dys_0"],Tuple [Scalar VDC DFloat "sm_0",Scalar VDC DFloat "u_0",Scalar VDC DFloat "v_0",Scalar VDC DFloat "w_0",Scalar VDC DFloat "g_2"],Scalar VDC DFloat "g_3"]),
        ("les_map_240",  [Tuple [Scalar VDC DFloat "dzn_0",Scalar VDC DFloat "dx1_0",Scalar VDC DFloat "dy1_0",Scalar VDC DFloat "dzs_0"],Tuple [Scalar VDC DFloat "sm_0",Scalar VDC DFloat "u_0",Scalar VDC DFloat "w_0",Scalar VDC DFloat "v_0",Scalar VDC DFloat "h_2"],Scalar VDC DFloat "h_3"]),
        ("press_map_138",  [Scalar VDC DFloat "pav_1",Scalar VDC DFloat "p0_1",Scalar VDC DFloat "p0_2"]),
        ("press_map_65",  [Tuple [Scalar VDC DFloat "dx1_0",Scalar VDC DFloat "dy1_0",Scalar VDC DFloat "dzn_0",Scalar VDC DFloat "dt_0"],Tuple [Scalar VDC DFloat "u_0",Scalar VDC DFloat "v_0",Scalar VDC DFloat "w_0",Scalar VDC DFloat "f_3",Scalar VDC DFloat "g_3",Scalar VDC DFloat "h_3",Scalar VDC DFloat "rhs_0"],Scalar VDC DFloat "rhs_1"]),
        ("press_map_87",  [Scalar VDC DFloat "rhsav_1",Scalar VDC DFloat "rhs_1",Scalar VDC DFloat "rhs_2"]),
        ("press_map_97",  [Tuple [Scalar VDC DFloat "dzs_0",Scalar VDC DFloat "dys_0",Scalar VDC DFloat "dxs_0",Scalar VDC DInt "nrd_0"],Tuple [Scalar VDC DFloat "p0_0",Scalar VDC DFloat "rhs_2",Scalar VDC DFloat "p1_0"],Tuple [Scalar VDC DFloat "p0_1",Scalar VDC DFloat "p1_1"]]),
        ("press_reduce_129",  [Tuple [Scalar VDC DFloat "dx1_0",Scalar VDC DFloat "dy1_0",Scalar VDC DFloat "dzn_0"],Tuple [Scalar VDC DFloat "pav_0",Scalar VDC DFloat "pco_0"],Scalar VDC DFloat "p0_1",Tuple [Scalar VDC DFloat "pav_1",Scalar VDC DFloat "pco_1"]]),
        ("press_reduce_78",  [Tuple [Scalar VDC DFloat "dx1_0",Scalar VDC DFloat "dy1_0",Scalar VDC DFloat "dzn_0"],Tuple [Scalar VDC DFloat "rhsav_0",Scalar VDC DFloat "area_0"],Scalar VDC DFloat "rhs_1",Tuple [Scalar VDC DFloat "area_1",Scalar VDC DFloat "rhsav_1"]]),
        ("velfg_map_135",  [Tuple [Scalar VDC DFloat "dy1_0",Scalar VDC DFloat "dx1_0",Scalar VDC DFloat "dzs_0",Scalar VDC DFloat "dzn_0"],Tuple [Scalar VDC DFloat "u_0",Scalar VDC DFloat "v_0",Scalar VDC DFloat "w_0"],Scalar VDC DFloat "g_4"]),
        ("velfg_map_175",  [Tuple [Scalar VDC DFloat "dzn_0",Scalar VDC DFloat "dx1_0",Scalar VDC DFloat "dy1_0"],Tuple [Scalar VDC DFloat "u_0",Scalar VDC DFloat "w_0",Scalar VDC DFloat "v_0"],Scalar VDC DFloat "h_4"]),
        ("velfg_map_95",  [Tuple [Scalar VDC DFloat "dx1_0",Scalar VDC DFloat "dy1_0",Scalar VDC DFloat "dzs_0",Scalar VDC DFloat "dzn_0"],Tuple [Scalar VDC DFloat "u_0",Scalar VDC DFloat "v_0",Scalar VDC DFloat "w_0"],Scalar VDC DFloat "f_4"]),
        ("velnw_map_36",  [Tuple [Scalar VDC DFloat "ro_0",Scalar VDC DFloat "dxs_0",Scalar VDC DFloat "dt_0"],Tuple [Scalar VDC DFloat "p0_2",Scalar VDC DFloat "u_0",Scalar VDC DFloat "f_4"],Scalar VDC DFloat "u_1"]),
        ("velnw_map_44",  [Tuple [Scalar VDC DFloat "ro_0",Scalar VDC DFloat "dys_0",Scalar VDC DFloat "dt_0"],Tuple [Scalar VDC DFloat "p0_2",Scalar VDC DFloat "v_0",Scalar VDC DFloat "g_4"],Scalar VDC DFloat "v_1"]),
        ("velnw_map_52",  [Tuple [Scalar VDC DFloat "ro_0",Scalar VDC DFloat "dzs_0",Scalar VDC DFloat "dt_0"],Tuple [Scalar VDC DFloat "p0_2",Scalar VDC DFloat "w_0",Scalar VDC DFloat "h_4"],Scalar VDC DFloat "w_1"])
    ]
stencilDefinitionsList = []

mainArgDeclsList = [
      ("bmask1_0" , MkFDecl "real"  (Just [7503492]) (Just In) ["bmask1_0"] )
    , ("fold_0" , MkFDecl "real"  (Just [7200000]) (Just In) ["fold_0"] )
    , ("wsum_0" , MkFDecl "real"  (Just [7338681]) (Just In) ["wsum_0"] )
    , ("dmask1_0" , MkFDecl "real"  (Just [7478728]) (Just In) ["dmask1_0"] )
    , ("dxs_0" , MkFDecl "real"  (Just [301]) (Just In) ["dxs_0"] )
    , ("pav_0" , MkFDecl "real" Nothing (Just In) ["pav_0"] )
    , ("w_0" , MkFDecl "real"  (Just [7594998]) (Just In) ["w_0"] )
    , ("dzs_0" , MkFDecl "real"  (Just [84]) (Just In) ["dzs_0"] )
    , ("dy1_0" , MkFDecl "real"  (Just [302]) (Just In) ["dy1_0"] )
    , ("p0_0" , MkFDecl "real"  (Just [7528338]) (Just In) ["p0_0"] )
    , ("f_0" , MkFDecl "real"  (Just [7338681]) (Just In) ["f_0"] )
    , ("dt_0" , MkFDecl "real" Nothing (Just In) ["dt_0"] )
    , ("cmask1_0" , MkFDecl "real"  (Just [7503492]) (Just In) ["cmask1_0"] )
    , ("dx1_0" , MkFDecl "real"  (Just [303]) (Just In) ["dx1_0"] )
    , ("usum_0" , MkFDecl "real"  (Just [7338681]) (Just In) ["usum_0"] )
    , ("beta_0" , MkFDecl "real" Nothing (Just In) ["beta_0"] )
    , ("v_0" , MkFDecl "real"  (Just [7503492]) (Just In) ["v_0"] )
    , ("g_0" , MkFDecl "real"  (Just [7338681]) (Just In) ["g_0"] )
    , ("hold_0" , MkFDecl "real"  (Just [7200000]) (Just In) ["hold_0"] )
    , ("h_0" , MkFDecl "real"  (Just [7338681]) (Just In) ["h_0"] )
    , ("delx1_0" , MkFDecl "real"  (Just [80]) (Just In) ["delx1_0"] )
    , ("alpha_0" , MkFDecl "real" Nothing (Just In) ["alpha_0"] )
    , ("dzn_0" , MkFDecl "real"  (Just [84]) (Just In) ["dzn_0"] )
    , ("gold_0" , MkFDecl "real"  (Just [7200000]) (Just In) ["gold_0"] )
    , ("dys_0" , MkFDecl "real"  (Just [301]) (Just In) ["dys_0"] )
    , ("rhsav_0" , MkFDecl "real" Nothing (Just In) ["rhsav_0"] )
    , ("nrd_0" , MkFDecl "integer" Nothing (Just In) ["nrd_0"] )
    , ("u_0" , MkFDecl "real"  (Just [7503492]) (Just In) ["u_0"] )
    , ("ro_0" , MkFDecl "real" Nothing (Just In) ["ro_0"] )
    , ("vsum_0" , MkFDecl "real"  (Just [7338681]) (Just In) ["vsum_0"] )
    , ("f_4" , MkFDecl "real"  (Just [7338681]) (Just Out) ["f_4"] )
    , ("sm_0" , MkFDecl "real"  (Just [7528338]) (Just Out) ["sm_0"] )
    , ("p1_1" , MkFDecl "real"  (Just [7528338]) (Just Out) ["p1_1"] )
    , ("fx_0" , MkFDecl "real"  (Just [7338681]) (Just Out) ["fx_0"] )
    , ("g_4" , MkFDecl "real"  (Just [7338681]) (Just Out) ["g_4"] )
    , ("fy_0" , MkFDecl "real"  (Just [7338681]) (Just Out) ["fy_0"] )
    , ("h_4" , MkFDecl "real"  (Just [7338681]) (Just Out) ["h_4"] )
    , ("fz_0" , MkFDecl "real"  (Just [7338681]) (Just Out) ["fz_0"] )
    , ("rhs_2" , MkFDecl "real"  (Just [7478728]) (Just Out) ["rhs_2"] )
  ]

scalarisedArgsList = [
     ( "adam_map_36",[("f",(0,InOut,"real")), ("g",(0,InOut,"real")), ("h",(0,InOut,"real")), ("fold",(0,InOut,"real")), ("gold",(0,InOut,"real")), ("hold",(0,InOut,"real"))])
    ,( "feedbf_map_49",[("usum",(0,InOut,"real")), ("u",(0,In,"real")), ("bmask1",(0,In,"real")), ("vsum",(0,InOut,"real")), ("v",(0,In,"real")), ("cmask1",(0,In,"real")), ("wsum",(0,InOut,"real")), ("w",(0,In,"real")), ("dmask1",(0,In,"real")), ("alpha",(0,In,"real")), ("dt",(0,In,"real")), ("beta",(0,In,"real")), ("fx",(0,Out,"real")), ("fy",(0,Out,"real")), ("fz",(0,Out,"real"))])
    ,( "feedbf_map_67",[("f",(0,InOut,"real")), ("fx",(0,In,"real")), ("g",(0,InOut,"real")), ("fy",(0,In,"real")), ("h",(0,InOut,"real")), ("fz",(0,In,"real"))])
    ,( "les_map_119",[("u",(0,In,"real")), ("dx1",(0,In,"real")), ("dy1",(0,In,"real")), ("dzs",(0,In,"real")), ("v",(0,In,"real")), ("w",(0,In,"real")), ("dzn",(0,In,"real")), ("delx1",(0,In,"real")), ("sm",(0,Out,"real"))])
    ,( "les_map_164",[("sm",(0,In,"real")), ("dy1",(0,In,"real")), ("dx1",(0,In,"real")), ("dzn",(0,In,"real")), ("u",(0,In,"real")), ("v",(0,In,"real")), ("dzs",(0,In,"real")), ("w",(0,In,"real")), ("dxs",(0,In,"real")), ("f",(0,InOut,"real"))])
    ,( "les_map_202",[("sm",(0,In,"real")), ("dy1",(0,In,"real")), ("dx1",(0,In,"real")), ("dzn",(0,In,"real")), ("u",(0,In,"real")), ("v",(0,In,"real")), ("dzs",(0,In,"real")), ("w",(0,In,"real")), ("dys",(0,In,"real")), ("g",(0,InOut,"real"))])
    ,( "les_map_240",[("sm",(0,In,"real")), ("dzn",(0,In,"real")), ("dx1",(0,In,"real")), ("dy1",(0,In,"real")), ("u",(0,In,"real")), ("dzs",(0,In,"real")), ("w",(0,In,"real")), ("v",(0,In,"real")), ("h",(0,InOut,"real"))])
    ,( "press_map_138",[("p0",(0,InOut,"real")), ("pav",(0,In,"real"))])
    ,( "press_map_65",[("u",(0,In,"real")), ("dx1",(0,In,"real")), ("v",(0,In,"real")), ("dy1",(0,In,"real")), ("w",(0,In,"real")), ("dzn",(0,In,"real")), ("f",(0,In,"real")), ("g",(0,In,"real")), ("h",(0,In,"real")), ("rhs",(0,InOut,"real")), ("dt",(0,In,"real"))])
    ,( "press_map_87",[("rhs",(0,InOut,"real")), ("rhsav",(0,In,"real"))])
    ,( "press_map_97",[("dzs",(0,In,"real")), ("dys",(0,In,"real")), ("dxs",(0,In,"real")), ("nrd",(0,In,"integer")), ("p0",(0,InOut,"real")), ("rhs",(0,In,"real")), ("p1",(0,InOut,"real"))])
    ,( "press_reduce_129",[("p0",(0,In,"real")), ("dx1",(0,In,"real")), ("dy1",(0,In,"real")), ("dzn",(0,In,"real")), ("pav",(0,InOut,"real")), ("pco",(0,InOut,"real"))])
    ,( "press_reduce_78",[("dx1",(0,In,"real")), ("dy1",(0,In,"real")), ("dzn",(0,In,"real")), ("rhs",(0,In,"real")), ("rhsav",(0,InOut,"real")), ("area",(0,InOut,"real"))])
    ,( "velfg_map_135",[("dy1",(0,In,"real")), ("u",(0,In,"real")), ("v",(0,In,"real")), ("dx1",(0,In,"real")), ("w",(0,In,"real")), ("dzs",(0,In,"real")), ("dzn",(0,In,"real")), ("g",(0,Out,"real"))])
    ,( "velfg_map_175",[("dzn",(0,In,"real")), ("u",(0,In,"real")), ("w",(0,In,"real")), ("dx1",(0,In,"real")), ("v",(0,In,"real")), ("dy1",(0,In,"real")), ("h",(0,Out,"real"))])
    ,( "velfg_map_95",[("u",(0,In,"real")), ("dx1",(0,In,"real")), ("v",(0,In,"real")), ("dy1",(0,In,"real")), ("w",(0,In,"real")), ("dzs",(0,In,"real")), ("dzn",(0,In,"real")), ("f",(0,Out,"real"))])
    ,( "velnw_map_36",[("p0",(0,In,"real")), ("ro",(0,In,"real")), ("dxs",(0,In,"real")), ("u",(0,In,"real")), ("dt",(0,In,"real")), ("f",(0,In,"real"))])
    ,( "velnw_map_44",[("p0",(0,In,"real")), ("ro",(0,In,"real")), ("dys",(0,In,"real")), ("v",(0,In,"real")), ("dt",(0,In,"real")), ("g",(0,In,"real"))])
    ,( "velnw_map_52",[("p0",(0,In,"real")), ("ro",(0,In,"real")), ("dzs",(0,In,"real")), ("w",(0,In,"real")), ("dt",(0,In,"real")), ("h",(0,In,"real"))])
  ]

origNamesList = [
     ( "adam_map_36",[("f",[("f_0",In),("f_1",Out)]), ("g",[("g_0",In),("g_1",Out)]), ("h",[("h_0",In),("h_1",Out)]), ("fold",[("fold_0",In),("fold_1",Out)]), ("gold",[("gold_0",In),("gold_1",Out)]), ("hold",[("hold_0",In),("hold_1",Out)])])
    ,( "feedbf_map_49",[("alpha",[("alpha_0",In)]), ("dt",[("dt_0",In)]), ("beta",[("beta_0",In)]), ("usum",[("usum_0",In),("usum_1",Out)]), ("u",[("u_0",In)]), ("bmask1",[("bmask1_0",In)]), ("vsum",[("vsum_0",In),("vsum_1",Out)]), ("v",[("v_0",In)]), ("cmask1",[("cmask1_0",In)]), ("wsum",[("wsum_0",In),("wsum_1",Out)]), ("w",[("w_0",In)]), ("dmask1",[("dmask1_0",In)]), ("fx",[("fx_0",Out)]), ("fy",[("fy_0",Out)]), ("fz",[("fz_0",Out)])])
    ,( "feedbf_map_67",[("f",[("f_1",In),("f_2",Out)]), ("fx",[("fx_0",In)]), ("g",[("g_1",In),("g_2",Out)]), ("fy",[("fy_0",In)]), ("h",[("h_1",In),("h_2",Out)]), ("fz",[("fz_0",In)])])
    ,( "les_map_119",[("dx1",[("dx1_0",In)]), ("dy1",[("dy1_0",In)]), ("dzs",[("dzs_0",In)]), ("dzn",[("dzn_0",In)]), ("delx1",[("delx1_0",In)]), ("u",[("u_0",In)]), ("v",[("v_0",In)]), ("w",[("w_0",In)]), ("sm",[("sm_0",Out)])])
    ,( "les_map_164",[("dy1",[("dy1_0",In)]), ("dx1",[("dx1_0",In)]), ("dzn",[("dzn_0",In)]), ("dzs",[("dzs_0",In)]), ("dxs",[("dxs_0",In)]), ("sm",[("sm_0",In)]), ("u",[("u_0",In)]), ("v",[("v_0",In)]), ("w",[("w_0",In)]), ("f",[("f_2",In),("f_3",Out)])])
    ,( "les_map_202",[("dy1",[("dy1_0",In)]), ("dx1",[("dx1_0",In)]), ("dzn",[("dzn_0",In)]), ("dzs",[("dzs_0",In)]), ("dys",[("dys_0",In)]), ("sm",[("sm_0",In)]), ("u",[("u_0",In)]), ("v",[("v_0",In)]), ("w",[("w_0",In)]), ("g",[("g_2",In),("g_3",Out)])])
    ,( "les_map_240",[("dzn",[("dzn_0",In)]), ("dx1",[("dx1_0",In)]), ("dy1",[("dy1_0",In)]), ("dzs",[("dzs_0",In)]), ("sm",[("sm_0",In)]), ("u",[("u_0",In)]), ("w",[("w_0",In)]), ("v",[("v_0",In)]), ("h",[("h_2",In),("h_3",Out)])])
    ,( "press_map_65",[("dx1",[("dx1_0",In)]), ("dy1",[("dy1_0",In)]), ("dzn",[("dzn_0",In)]), ("dt",[("dt_0",In)]), ("u",[("u_0",In)]), ("v",[("v_0",In)]), ("w",[("w_0",In)]), ("f",[("f_3",In)]), ("g",[("g_3",In)]), ("h",[("h_3",In)]), ("rhs",[("rhs_0",In),("rhs_1",Out)])])
    ,( "press_reduce_78",[("area",[("area_1",Out)]), ("rhsav",[("rhsav_1",Out)]), ("rhs",[("rhs_1",In)]), ("dx1",[("dx1_0",In)]), ("dy1",[("dy1_0",In)]), ("dzn",[("dzn_0",In)]), ("rhsav",[("rhsav_0",Out)]), ("area",[("area_0",Out)])])
    ,( "press_map_87",[("rhsav",[("rhsav_1",In)]), ("rhs",[("rhs_1",In),("rhs_2",Out)])])
    ,( "press_map_97",[("dzs",[("dzs_0",In)]), ("dys",[("dys_0",In)]), ("dxs",[("dxs_0",In)]), ("nrd",[("nrd_0",In)]), ("p0",[("p0_0",In),("p0_1",Out)]), ("rhs",[("rhs_2",In)]), ("p1",[("p1_0",In),("p1_1",Out)])])
    ,( "press_reduce_129",[("pav",[("pav_1",Out)]), ("pco",[("pco_1",Out)]), ("p0",[("p0_1",In)]), ("dx1",[("dx1_0",In)]), ("dy1",[("dy1_0",In)]), ("dzn",[("dzn_0",In)]), ("pav",[("pav_0",Out)]), ("pco",[("pco_0",Out)])])
    ,( "press_map_138",[("pav",[("pav_1",In)]), ("p0",[("p0_1",In),("p0_2",Out)])])
    ,( "velfg_map_95",[("dx1",[("dx1_0",In)]), ("dy1",[("dy1_0",In)]), ("dzs",[("dzs_0",In)]), ("dzn",[("dzn_0",In)]), ("u",[("u_0",In)]), ("v",[("v_0",In)]), ("w",[("w_0",In)]), ("f",[("f_4",Out)])])
    ,( "velfg_map_135",[("dy1",[("dy1_0",In)]), ("dx1",[("dx1_0",In)]), ("dzs",[("dzs_0",In)]), ("dzn",[("dzn_0",In)]), ("u",[("u_0",In)]), ("v",[("v_0",In)]), ("w",[("w_0",In)]), ("g",[("g_4",Out)])])
    ,( "velfg_map_175",[("dzn",[("dzn_0",In)]), ("dx1",[("dx1_0",In)]), ("dy1",[("dy1_0",In)]), ("u",[("u_0",In)]), ("w",[("w_0",In)]), ("v",[("v_0",In)]), ("h",[("h_4",Out)])])
    ,( "velnw_map_36",[("ro",[("ro_0",In)]), ("dxs",[("dxs_0",In)]), ("dt",[("dt_0",In)]), ("p0",[("p0_2",In)]), ("u",[("u_0",In),("u_1",Out)]), ("f",[("f_4",In)])])
    ,( "velnw_map_44",[("ro",[("ro_0",In)]), ("dys",[("dys_0",In)]), ("dt",[("dt_0",In)]), ("p0",[("p0_2",In)]), ("v",[("v_0",In),("v_1",Out)]), ("g",[("g_4",In)])])
    ,( "velnw_map_52",[("ro",[("ro_0",In)]), ("dzs",[("dzs_0",In)]), ("dt",[("dt_0",In)]), ("p0",[("p0_2",In)]), ("w",[("w_0",In),("w_1",Out)]), ("h",[("h_4",In)])])
  ]

{-
f_0 :: Vec 7338681 Float
g_0 :: Vec 7338681 Float
h_0 :: Vec 7338681 Float
fold_0 :: Vec 7200000 Float
gold_0 :: Vec 7200000 Float
hold_0 :: Vec 7200000 Float
usum_0 :: Vec 7338681 Float
u_0 :: Vec 7503492 Float
bmask1_0 :: Vec 7503492 Float
vsum_0 :: Vec 7338681 Float
v_0 :: Vec 7503492 Float
cmask1_0 :: Vec 7503492 Float
wsum_0 :: Vec 7338681 Float
w_0 :: Vec 7594998 Float
dmask1_0 :: Vec 7478728 Float
alpha_0 :: Float
dt_0 :: Float
beta_0 :: Float
dx1_0 :: SVec 303 Float
dy1_0 :: SVec 302 Float
dzs_0 :: SVec 84 Float
dzn_0 :: SVec 84 Float
delx1_0 :: SVec 80 Float
dxs_0 :: SVec 301 Float
dys_0 :: SVec 301 Float
rhsav_0 :: Float
p0_0 :: Vec 7528338 Float
nrd_0 :: Int
pav_0 :: Float
ro_0 :: Float

fx_0 :: Vec 7338681 Float
fy_0 :: Vec 7338681 Float
fz_0 :: Vec 7338681 Float
sm_0 :: Vec 7528338 Float
rhs_2 :: Vec 7478728 Float
p1_1 :: Vec 7528338 Float
f_4 :: Vec 7338681 Float
g_4 :: Vec 7338681 Float
h_4 :: Vec 7338681 Float

velnw_map_52 :: (Float,Float,Float) -> (Float,Float,Float) -> Float
velnw_map_44 :: (Float,Float,Float) -> (Float,Float,Float) -> Float
velnw_map_36 :: (Float,Float,Float) -> (Float,Float,Float) -> Float
velfg_map_95 :: (Float,Float,Float,Float) -> (Float,Float,Float) -> Float
velfg_map_175 :: (Float,Float,Float) -> (Float,Float,Float) -> Float
velfg_map_135 :: (Float,Float,Float,Float) -> (Float,Float,Float) -> Float
press_reduce_78 :: (Float,Float,Float) -> (Float,Float) -> Float -> (Float,Float)
press_reduce_129 :: (Float,Float,Float) -> (Float,Float) -> Float -> (Float,Float)
press_map_97 :: (Float,Float,Float,Int) -> (Float,Float,Float) -> (Float,Float)
press_map_87 :: Float -> Float -> Float
press_map_65 :: (Float,Float,Float,Float) -> (Float,Float,Float,Float,Float,Float,Float) -> Float
press_map_138 :: Float -> Float -> Float
les_map_240 :: (Float,Float,Float,Float) -> (Float,Float,Float,Float,Float) -> Float
les_map_202 :: (Float,Float,Float,Float,Float) -> (Float,Float,Float,Float,Float) -> Float
les_map_164 :: (Float,Float,Float,Float,Float) -> (Float,Float,Float,Float,Float) -> Float
les_map_119 :: (Float,Float,Float,Float,Float) -> (Float,Float,Float) -> Float
feedbf_map_67 :: (Float,Float,Float,Float,Float,Float) -> (Float,Float,Float)
feedbf_map_49 :: (Float,Float,Float) -> (Float,Float,Float,Float,Float,Float,Float,Float,Float) -> (Float,Float,Float,Float,Float,Float)
adam_map_36 :: (Float,Float,Float,Float,Float,Float) -> (Float,Float,Float,Float,Float,Float)

main :: (Vec 7338681 Float,Vec 7338681 Float,Vec 7338681 Float,Vec 7200000 Float,Vec 7200000 Float,Vec 7200000 Float,Vec 7338681 Float,Vec 7503492 Float,Vec 7503492 Float,Vec 7338681 Float,Vec 7503492 Float,Vec 7503492 Float,Vec 7338681 Float,Vec 7594998 Float,Vec 7478728 Float,Float,Float,Float,SVec 303 Float,SVec 302 Float,SVec 84 Float,SVec 84 Float,SVec 80 Float,SVec 301 Float,SVec 301 Float,Float,Vec 7528338 Float,Int,Float,Float) -> (Vec 7338681 Float,Vec 7338681 Float,Vec 7338681 Float,Vec 7528338 Float,Vec 7478728 Float,Vec 7528338 Float,Vec 7338681 Float,Vec 7338681 Float,Vec 7338681 Float)
main (f_0,g_0,h_0,fold_0,gold_0,hold_0,usum_0,u_0,bmask1_0,vsum_0,v_0,cmask1_0,wsum_0,w_0,dmask1_0,alpha_0,dt_0,beta_0,dx1_0,dy1_0,dzs_0,dzn_0,delx1_0,dxs_0,dys_0,rhsav_0,p0_0,nrd_0,pav_0,ro_0) =
  let
     -- adam_map_36
    (f_1,fold_1,g_1,gold_1,h_1,hold_1) = unzipt $ map adam_map_36 (zipt (f_0,g_0,h_0,fold_0,gold_0,hold_0))
     -- feedbf_map_49
    (fx_0,fy_0,fz_0,usum_1,vsum_1,wsum_1) = unzipt $ map (feedbf_map_49 (alpha_0,dt_0,beta_0)) (zipt (usum_0,u_0,bmask1_0,vsum_0,v_0,cmask1_0,wsum_0,w_0,dmask1_0))
     -- feedbf_map_67
    (f_2,g_2,h_2) = unzipt $ map feedbf_map_67 (zipt (f_1,fx_0,g_1,fy_0,h_1,fz_0))
     -- les_map_119
    sm_0 =  map (les_map_119 (dx1_0,dy1_0,dzs_0,dzn_0,delx1_0)) (zipt (u_0,v_0,w_0))
     -- les_map_164
    f_3 =  map (les_map_164 (dy1_0,dx1_0,dzn_0,dzs_0,dxs_0)) (zipt (sm_0,u_0,v_0,w_0,f_2))
     -- les_map_202
    g_3 =  map (les_map_202 (dy1_0,dx1_0,dzn_0,dzs_0,dys_0)) (zipt (sm_0,u_0,v_0,w_0,g_2))
     -- les_map_240
    h_3 =  map (les_map_240 (dzn_0,dx1_0,dy1_0,dzs_0)) (zipt (sm_0,u_0,w_0,v_0,h_2))
     -- press_map_65
    rhs_1 =  map (press_map_65 (dx1_0,dy1_0,dzn_0,dt_0)) (zipt (u_0,v_0,w_0,f_3,g_3,h_3,rhs_0))
     -- press_reduce_78
    (area_1,rhsav_1) =  fold (press_reduce_78 (dx1_0,dy1_0,dzn_0)) (rhsav_0,area_0) rhs_1
     -- press_map_87
    rhs_2 =  map (press_map_87 rhsav_1) rhs_1
     -- press_map_97
    (p0_1,p1_1) = unzipt $ map (press_map_97 (dzs_0,dys_0,dxs_0,nrd_0)) (zipt (p0_0,rhs_2,p1_0))
     -- press_reduce_129
    (pav_1,pco_1) =  fold (press_reduce_129 (dx1_0,dy1_0,dzn_0)) (pav_0,pco_0) p0_1
     -- press_map_138
    p0_2 =  map (press_map_138 pav_1) p0_1
     -- velfg_map_95
    f_4 =  map (velfg_map_95 (dx1_0,dy1_0,dzs_0,dzn_0)) (zipt (u_0,v_0,w_0))
     -- velfg_map_135
    g_4 =  map (velfg_map_135 (dy1_0,dx1_0,dzs_0,dzn_0)) (zipt (u_0,v_0,w_0))
     -- velfg_map_175
    h_4 =  map (velfg_map_175 (dzn_0,dx1_0,dy1_0)) (zipt (u_0,w_0,v_0))
     -- velnw_map_36
    u_1 =  map (velnw_map_36 (ro_0,dxs_0,dt_0)) (zipt (p0_2,u_0,f_4))
     -- velnw_map_44
    v_1 =  map (velnw_map_44 (ro_0,dys_0,dt_0)) (zipt (p0_2,v_0,g_4))
     -- velnw_map_52
    w_1 =  map (velnw_map_52 (ro_0,dzs_0,dt_0)) (zipt (p0_2,w_0,h_4))
  in
    (fx_0,fy_0,fz_0,sm_0,rhs_2,p1_1,f_4,g_4,h_4)
-}
