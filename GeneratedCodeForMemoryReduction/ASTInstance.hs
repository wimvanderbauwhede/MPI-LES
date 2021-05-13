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
    , ("g_0" , MkFDecl "real"  (Just [7338681]) (Just In) ["g_0"] )
    , ("dt_0" , MkFDecl "real" Nothing (Just In) ["dt_0"] )
    , ("f_0" , MkFDecl "real"  (Just [7338681]) (Just In) ["f_0"] )
    , ("pav_0" , MkFDecl "real" Nothing (Just In) ["pav_0"] )
    , ("nrd_0" , MkFDecl "integer" Nothing (Just In) ["nrd_0"] )
    , ("dzn_0" , MkFDecl "real"  (Just [84]) (Just In) ["dzn_0"] )
    , ("usum_0" , MkFDecl "real"  (Just [7338681]) (Just In) ["usum_0"] )
    , ("w_0" , MkFDecl "real"  (Just [7594998]) (Just In) ["w_0"] )
    , ("p0_0" , MkFDecl "real"  (Just [7528338]) (Just In) ["p0_0"] )
    , ("gold_0" , MkFDecl "real"  (Just [7200000]) (Just In) ["gold_0"] )
    , ("dy1_0" , MkFDecl "real"  (Just [302]) (Just In) ["dy1_0"] )
    , ("dxs_0" , MkFDecl "real"  (Just [301]) (Just In) ["dxs_0"] )
    , ("fold_0" , MkFDecl "real"  (Just [7200000]) (Just In) ["fold_0"] )
    , ("v_0" , MkFDecl "real"  (Just [7503492]) (Just In) ["v_0"] )
    , ("h_0" , MkFDecl "real"  (Just [7338681]) (Just In) ["h_0"] )
    , ("vsum_0" , MkFDecl "real"  (Just [7338681]) (Just In) ["vsum_0"] )
    , ("ro_0" , MkFDecl "real" Nothing (Just In) ["ro_0"] )
    , ("u_0" , MkFDecl "real"  (Just [7503492]) (Just In) ["u_0"] )
    , ("dzs_0" , MkFDecl "real"  (Just [84]) (Just In) ["dzs_0"] )
    , ("wsum_0" , MkFDecl "real"  (Just [7338681]) (Just In) ["wsum_0"] )
    , ("hold_0" , MkFDecl "real"  (Just [7200000]) (Just In) ["hold_0"] )
    , ("delx1_0" , MkFDecl "real"  (Just [80]) (Just In) ["delx1_0"] )
    , ("beta_0" , MkFDecl "real" Nothing (Just In) ["beta_0"] )
    , ("cmask1_0" , MkFDecl "real"  (Just [7503492]) (Just In) ["cmask1_0"] )
    , ("dmask1_0" , MkFDecl "real"  (Just [7478728]) (Just In) ["dmask1_0"] )
    , ("dx1_0" , MkFDecl "real"  (Just [303]) (Just In) ["dx1_0"] )
    , ("dys_0" , MkFDecl "real"  (Just [301]) (Just In) ["dys_0"] )
    , ("rhsav_0" , MkFDecl "real" Nothing (Just In) ["rhsav_0"] )
    , ("alpha_0" , MkFDecl "real" Nothing (Just In) ["alpha_0"] )
    , ("fz_0" , MkFDecl "real"  (Just [7338681]) (Just Out) ["fz_0"] )
    , ("f_4" , MkFDecl "real"  (Just [7338681]) (Just Out) ["f_4"] )
    , ("p1_1" , MkFDecl "real"  (Just [7528338]) (Just Out) ["p1_1"] )
    , ("sm_0" , MkFDecl "real"  (Just [7528338]) (Just Out) ["sm_0"] )
    , ("fy_0" , MkFDecl "real"  (Just [7338681]) (Just Out) ["fy_0"] )
    , ("g_4" , MkFDecl "real"  (Just [7338681]) (Just Out) ["g_4"] )
    , ("rhs_2" , MkFDecl "real"  (Just [7478728]) (Just Out) ["rhs_2"] )
    , ("h_4" , MkFDecl "real"  (Just [7338681]) (Just Out) ["h_4"] )
    , ("fx_0" , MkFDecl "real"  (Just [7338681]) (Just Out) ["fx_0"] )
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