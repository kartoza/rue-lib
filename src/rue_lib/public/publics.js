/**
 * to use this code: import publics from this js file as well as the GI module
 * run publics with the GI module as input along with other start node input
 * e.g.:
 * const publics = require('./publics.js').publics
 * const module = require('gi-module')
 * const result = await publics(module, start_input_1, start_input_2, ...);
 *
 * returns: a json object:
 *   _ result.model -> gi model of the flowchart
 *   _ result.result -> returned output of the flowchart, if the flowchart does not return any value, result.result is the model of the flowchart
 */

// Parameter: {"name":"IN_MODEL","value":"cluster.sim","type":0}
// Parameter: {"name":"FORCE_STOP_MOBIUS","value":"[false]","type":0}
// Parameter: {"name":"ROAD_ART_W","value":"20","type":1,"min":"10","max":"30","step":"1"}
// Parameter: {"name":"ROAD_SEC_W","value":"15","type":1,"min":"10","max":"30","step":"1"}
// Parameter: {"name":"ROAD_LOC_W","value":10,"type":1,"min":"5","max":"30","step":"1"}
// Parameter: {"name":"PART_ART_D","value":"40","type":1,"min":"5","max":"50","step":"1"}
// Parameter: {"name":"PART_SEC_D","value":30,"type":1,"min":"5","max":"50","step":"1"}
// Parameter: {"name":"PART_LOC_D","value":20,"type":1,"min":"5","max":"50","step":"1"}
// Parameter: {"name":"PART_OG_D","value":"40","type":1,"min":"5","max":"50","step":"1"}
// Parameter: {"name":"PART_OG_W","value":30,"type":1,"min":"5","max":"50","step":"1"}
// Parameter: {"name":"PLOT_ART_W","value":40,"type":1,"min":"5","max":"50","step":"1"}
// Parameter: {"name":"PLOT_SEC_W","value":20,"type":1,"min":"5","max":"50","step":"3"}
// Parameter: {"name":"PLOT_LOC_W","value":11,"type":1,"min":"5","max":"50","step":"3"}
// Parameter: {"name":"BLK_ART_NUM_OG_D","value":0,"type":1,"min":"0","max":"1","step":"1"}
// Parameter: {"name":"BLK_ART_NUM_OG_W","value":3,"type":1,"min":"0","max":"10","step":"1"}
// Parameter: {"name":"BLK_SEC_NUM_OG_D","value":0,"type":1,"min":"0","max":"1","step":"1"}
// Parameter: {"name":"BLK_SEC_NUM_OG_W","value":3,"type":1,"min":"0","max":"10","step":"1"}
// Parameter: {"name":"BLK_LOC_NUM_OG_D","value":2,"type":1,"min":"0","max":"2","step":"1"}
// Parameter: {"name":"BLK_LOC_NUM_OG_W","value":3,"type":1,"min":"0","max":"10","step":"1"}
// Parameter: {"name":"PATH_W","value":4,"type":1,"min":"3","max":"10","step":"1"}
// Parameter: {"name":"OPEN_PERCENT","value":"6","type":1,"min":"0","max":"10","step":"0.1"}
// Parameter: {"name":"AMEN_PERCENT","value":"8","type":1,"min":"0","max":"15","step":"0.1"}
// Parameter: {"name":"PAVEMENT_W","value":3,"type":1,"min":"1","max":"10","step":"0.5"}
// Parameter: {"name":"ADD_TREES","value":true,"type":2}
// Parameter: {"name":"TREE_SPACING","value":12,"type":1,"min":"6","max":"30","step":"1"}
// Parameter: {"name":"TREE_HEIGHT_START","value":"8","type":1,"min":"5","max":"15","step":"1"}
// Parameter: {"name":"TREE_HEIGHT_MAX","value":"20","type":1,"min":"10","max":"30","step":"1"}


const mfn = require('@design-automation/mobius-sim-funcs').Funcs();
const ifn = require('@design-automation/mobius-inline-funcs').InlineClass(true);

async function publics(IN_MODEL, FORCE_STOP_MOBIUS, ROAD_ART_W, ROAD_SEC_W, ROAD_LOC_W, PART_ART_D, PART_SEC_D, PART_LOC_D, PART_OG_D, PART_OG_W, PLOT_ART_W, PLOT_SEC_W, PLOT_LOC_W, BLK_ART_NUM_OG_D, BLK_ART_NUM_OG_W, BLK_SEC_NUM_OG_D, BLK_SEC_NUM_OG_W, BLK_LOC_NUM_OG_D, BLK_LOC_NUM_OG_W, PATH_W, OPEN_PERCENT, AMEN_PERCENT, PAVEMENT_W, ADD_TREES, TREE_SPACING, TREE_HEIGHT_START, TREE_HEIGHT_MAX) {

  var __model__ = null;

  /** * **/

  async function exec_publics_Zone($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {


    async function exec_publics_Zone_node_xlp2w5q5fr_rotSite_($p, pline_art_site_, site_) {

      let posis_ = mfn.query.Get('ps', pline_art_site_);

      let straight_pline_ = mfn.make.Polyline([posis_[pythonList(0, posis_.length)], posis_[pythonList(-1, posis_.length)]], 'open');

      let ray_ = mfn.calc.Ray(mfn.query.Get('_e', straight_pline_)[pythonList(0, mfn.query.Get('_e', straight_pline_).length)]);

      let vecx_ = ifn.vecNorm(ray_[pythonList(1, ray_.length)]);

      let vecy_ = [-vecx_[pythonList(1, vecx_.length)], vecx_[pythonList(0, vecx_.length)], 0];

      let pln_ = ifn.plnMake(ray_[pythonList(0, ray_.length)], vecx_, vecy_);

      mfn.modify.XForm(site_, pln_, JSON.parse(JSON.stringify(ifn.XY)));

      mfn.edit.Delete(straight_pline_, 'delete_selected');

      return pln_;
    }


    async function exec_publics_Zone_node_xlp2w5q5fr_getPerimPlines_($p, site_, road_descr_) {

      let edges_ = [];

      for (let edge_ of mfn.query.Get('_e', site_)) {

        if (road_descr_ == mfn.attrib.Get(edge_, 'road')) {

          mfn.list.Add(edges_, edge_, 'to_end');
        }
      }

      let posis_ = [mfn.query.Get('ps', edges_[pythonList(0, edges_.length)])];

      for (let edge_ of edges_.slice(1)) {

        let start_end_ = mfn.query.Get('ps', edge_);

        if (start_end_[pythonList(0, start_end_.length)] == posis_[pythonList(-1, posis_.length)][pythonList(-1, posis_[pythonList(-1, posis_.length)].length)]) {

          mfn.list.Add(posis_[pythonList(-1, posis_.length)], start_end_[pythonList(1, start_end_.length)], 'to_end');
        } else {
          if (start_end_[pythonList(1, start_end_.length)] == posis_[pythonList(-1, posis_.length)][pythonList(0, posis_[pythonList(-1, posis_.length)].length)]) {

            mfn.list.Add(posis_[pythonList(-1, posis_.length)], start_end_[pythonList(0, start_end_.length)], 'to_start');
          } else {

            mfn.list.Add(posis_, start_end_, 'to_end');
          }
        }
      }

      let site_plines_ = mfn.make.Polyline(posis_, 'open');

      return site_plines_;
    }


    async function exec_publics_Zone_node_xlp2w5q5fr_extendPline_($p, plines_, dist_) {

      for (let pline_ of plines_) {

        let closed_ = mfn.query.Type(pline_, 'is_closed');

        if (!closed_) {

          let var_ = mfn.edit.Weld(pline_, 'break_weld');

          let edges_ = mfn.query.Get('_e', pline_);

          let posis_ = mfn.query.Get('ps', pline_);

          let xyzs_ = mfn.attrib.Get(mfn.query.Get('ps', edges_[pythonList(0, edges_.length)]), 'xyz');

          let vec_ = ifn.vecSetLen(ifn.vecFromTo(xyzs_[pythonList(1, xyzs_.length)], xyzs_[pythonList(0, xyzs_.length)]), dist_);

          mfn.modify.Move(posis_[pythonList(0, posis_.length)], vec_);

          xyzs_ = mfn.attrib.Get(mfn.query.Get('ps', edges_[pythonList(-1, edges_.length)]), 'xyz');

          vec_ = ifn.vecSetLen(ifn.vecFromTo(xyzs_[pythonList(0, xyzs_.length)], xyzs_[pythonList(1, xyzs_.length)]), dist_);

          mfn.modify.Move(posis_[pythonList(-1, posis_.length)], vec_);
        }
      }
    }

    async function exec_publics_Zone_node_xlp2w5q5fr($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: rot', '__null__')
      }


      let coords_ = SITE_[pythonList(0, SITE_.length)];

      let edgesroad_ = SITE_[pythonList(1, SITE_.length)];

      let ps_ = mfn.make.Position(coords_);

      let pg_ = mfn.make.Polygon(ps_);

      let site1_ = pg_;

      let edges_ = mfn.query.Get('_e', site1_);

      let p_ = mfn.attrib.Get(mfn.query.Get('ps', site1_), 'xyz');

      for (let i_ of ifn.range(ifn.len(edges_))) {

        let e_ = edges_[pythonList(i_, edges_.length)];

        mfn.attrib.Set(e_, `road`, edgesroad_[pythonList(i_, edgesroad_.length)]);
      }

      mfn.attrib.Set(site1_, `type`, 'main_site');

      let edges0_ = mfn.query.Get('_e', site1_);

      let edges1_ = mfn.query.Get('_e', site1_);

      for (let i_ of ifn.range(ifn.len(edges0_))) {

        mfn.attrib.Set(edges1_[pythonList(i_, edges1_.length)], `road`, mfn.attrib.Get(edges0_[pythonList(i_, edges0_.length)], 'road'));
      }

      let pline_art_site_ = await exec_publics_Zone_node_xlp2w5q5fr_getPerimPlines_($p, site1_, "road_art");
      if ($p.terminated) {
        return mfn.getModel();
      }

      mfn.attrib.Set(null, `pln`, await exec_publics_Zone_node_xlp2w5q5fr_rotSite_($p, pline_art_site_, site1_));
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_publics_Zone_node_xlp2w5q5fr_extendPline_($p, pline_art_site_, 200);
      if ($p.terminated) {
        return mfn.getModel();
      }

      mfn.attrib.Set(pline_art_site_, `name`, "art_pline");
    }


    async function exec_publics_Zone_node_9ylsi2n679i_getLongestEdge_($p, edge_) {

      let road_ = mfn.attrib.Get(edge_, 'road');

      let edges_same_type_ = [edge_];

      let this_edge_ = edge_;

      let this_road_ = road_;

      while (this_road_ == road_) {

        let force_stop_mobius0_ = 'check_force_stop';

        this_edge_ = await exec_publics_Zone_node_9ylsi2n679i_nextEdge_($p, this_edge_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        this_road_ = mfn.attrib.Get(this_edge_, 'road');

        if (this_road_ == road_) {

          mfn.list.Add(edges_same_type_, this_edge_, 'to_end');
        }
      }

      this_edge_ = edge_;

      this_road_ = road_;

      while (this_road_ == road_) {

        let force_stop_mobius0_ = 'check_force_stop';

        this_edge_ = await exec_publics_Zone_node_9ylsi2n679i_prevEdge_($p, this_edge_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        this_road_ = mfn.attrib.Get(this_edge_, 'road');

        if (this_road_ == road_) {

          mfn.list.Add(edges_same_type_, this_edge_, 'to_end');
        }
      }

      edges_same_type_ = edges_same_type_;

      let lengths_ = mfn.calc.Length(edges_same_type_);

      return ifn.listSort(edges_same_type_, lengths_)[pythonList(-1, ifn.listSort(edges_same_type_, lengths_).length)];
    }


    async function exec_publics_Zone_node_9ylsi2n679i_getSiteLongestEdges_($p, site_, road_descr_) {

      let edges_ = [];

      for (let edge_ of mfn.query.Get('_e', site_)) {

        if (road_descr_ == mfn.attrib.Get(edge_, 'road')) {

          mfn.list.Add(edges_, edge_, 'to_end');
        }
      }

      if (ifn.len(edges_) == 0) {

        return [];
      }

      let posis_ = [mfn.query.Get('ps', edges_[pythonList(0, edges_.length)])];

      for (let edge_ of edges_.slice(1)) {

        let start_end_ = mfn.query.Get('ps', edge_);

        if (start_end_[pythonList(0, start_end_.length)] == posis_[pythonList(-1, posis_.length)][pythonList(-1, posis_[pythonList(-1, posis_.length)].length)]) {

          mfn.list.Add(posis_[pythonList(-1, posis_.length)], start_end_[pythonList(1, start_end_.length)], 'to_end');
        } else {
          if (start_end_[pythonList(1, start_end_.length)] == posis_[pythonList(-1, posis_.length)][pythonList(0, posis_[pythonList(-1, posis_.length)].length)]) {

            mfn.list.Add(posis_[pythonList(-1, posis_.length)], start_end_[pythonList(0, start_end_.length)], 'to_start');
          } else {

            mfn.list.Add(posis_, start_end_, 'to_end');
          }
        }
      }

      let site_plines_ = mfn.make.Polyline(posis_, 'open');

      return site_plines_;
    }


    async function exec_publics_Zone_node_9ylsi2n679i_extend_($p, pline_, dist_) {

      let edges_ = mfn.query.Get('_e', pline_);

      let lengths_ = mfn.calc.Length(edges_);

      let longest_ = ifn.listSort(edges_, lengths_)[pythonList(-1, ifn.listSort(edges_, lengths_).length)];

      let xyzs_ = mfn.attrib.Get(mfn.query.Get('ps', longest_), 'xyz');

      let vec_ = ifn.vecSetLen(ifn.vecFromTo(xyzs_[pythonList(0, xyzs_.length)], xyzs_[pythonList(1, xyzs_.length)]), dist_);

      let xyz0_ = ifn.vecSub(xyzs_[pythonList(0, xyzs_.length)], vec_);

      let xyz1_ = ifn.vecAdd(xyzs_[pythonList(1, xyzs_.length)], vec_);

      let posis_ = mfn.make.Position([xyz0_, xyz1_]);

      let long_pline_ = mfn.make.Polyline(posis_, 'open');

      return long_pline_;
    }


    async function exec_publics_Zone_node_9ylsi2n679i_offsetEdge_($p, site_, pline_, dist_) {

      let pline_long_ = await exec_publics_Zone_node_9ylsi2n679i_extend_($p, pline_, 1000);
      if ($p.terminated) {
        return mfn.getModel();
      }

      mfn.modify.Offset(pline_long_, -dist_);

      pline_ = await exec_publics_Zone_node_9ylsi2n679i_siteTrimPline_($p, site_, pline_long_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      return pline_;
    }


    async function exec_publics_Zone_node_9ylsi2n679i_siteTrimPline_($p, site_, pline_) {

      let pline2_ = mfn.poly2d.Boolean(pline_, site_, 'intersect');

      mfn.edit.Delete(pline_, 'delete_selected');

      return pline2_;
    }


    async function exec_publics_Zone_node_9ylsi2n679i_getEdge_($p, wire_, posi_) {

      let posis_ = mfn.query.Get('ps', wire_);

      let edges_ = mfn.query.Get('_e', wire_);

      let idx_ = ifn.listFind(posis_, posi_);

      if (idx_ == -1) {

        return null;
      }

      return edges_[pythonList(idx_, edges_.length)];
    }


    async function exec_publics_Zone_node_9ylsi2n679i_nextEdge_($p, edge_) {

      let wire_ = mfn.query.Get('_w', edge_);

      let edges_ = mfn.query.Get('_e', wire_);

      let idx_ = ifn.listFind(edges_, edge_);

      let next_idx_ = idx_ + 1;

      if (next_idx_ == ifn.len(edges_)) {

        return edges_[pythonList(0, edges_.length)];
      }

      return edges_[pythonList(next_idx_, edges_.length)];
    }


    async function exec_publics_Zone_node_9ylsi2n679i_prevEdge_($p, edge_) {

      let wire_ = mfn.query.Get('_w', edge_);

      let edges_ = mfn.query.Get('_e', wire_);

      let idx_ = ifn.listFind(edges_, edge_);

      let prev_idx_ = idx_ - 1;

      if (prev_idx_ < 0) {

        return edges_[pythonList(-1, edges_.length)];
      }

      return edges_[pythonList(prev_idx_, edges_.length)];
    }


    async function exec_publics_Zone_node_9ylsi2n679i_sortAlongX_($p, plines_) {

      let cens_ = mfn.calc.Centroid(plines_, 'ps_average');

      let x_vals_ = [];

      for (let cen_ of cens_) {

        mfn.list.Add(x_vals_, cen_[pythonList(0, cen_.length)], 'to_end');
      }

      let sorted_ = ifn.listSort(ifn.listZip([x_vals_, plines_]), 0);

      let sorted_plines_ = [];

      for (let x_pline_ of sorted_) {

        mfn.list.Add(sorted_plines_, x_pline_[pythonList(1, x_pline_.length)], 'to_end');
      }

      return sorted_plines_;
    }


    async function exec_publics_Zone_node_9ylsi2n679i_closestEdge_($p, site_, xyz_) {

      let min_d_ = Infinity;

      let closest_edge_ = null;

      if (xyz_[pythonList(0, xyz_.length)] > 1290) {
      }

      for (let edge_ of mfn.query.Get('_e', site_)) {

        let r_ = mfn.calc.Ray(edge_);

        let d_ = ifn.distance(xyz_, ifn.project(xyz_, r_, 0));

        if (d_ < min_d_) {

          min_d_ = d_;

          closest_edge_ = edge_;
        }

        if (d_ < 0.001) {

          break;
        }
      }

      return closest_edge_;
    }

    async function exec_publics_Zone_node_9ylsi2n679i($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: left and right plines', '__null__')
      }


      let site1_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', 'main_site')[pythonList(0, 'main_site'.length)];

      let site_art_ = (mfn.query.Filter(mfn.query.Get('_e', null), ['road', null], '==', "road_art"));

      let posis_ = mfn.query.Get('ps', site_art_);

      let edge0_ = await exec_publics_Zone_node_9ylsi2n679i_getEdge_($p, site1_, posis_[pythonList(0, posis_.length)]);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let left_edge_ = await exec_publics_Zone_node_9ylsi2n679i_prevEdge_($p, edge0_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let left_type_ = mfn.attrib.Get(left_edge_, 'road');

      if (left_type_ != "cold") {

        let left_edge_long_ = await exec_publics_Zone_node_9ylsi2n679i_getLongestEdge_($p, left_edge_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        let left_pline_ = await exec_publics_Zone_node_9ylsi2n679i_extend_($p, left_edge_long_, 500);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(left_pline_, `name`, "left_pline");

        mfn.attrib.Set(left_pline_, `type`, left_type_);

        mfn.attrib.Set(left_pline_, `site_edge`, left_edge_long_);
      }

      let right_edge_ = await exec_publics_Zone_node_9ylsi2n679i_getEdge_($p, site1_, posis_[pythonList(-1, posis_.length)]);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let right_type_ = mfn.attrib.Get(right_edge_, 'road');

      if (right_type_ != "cold") {

        let right_edge_long_ = await exec_publics_Zone_node_9ylsi2n679i_getLongestEdge_($p, right_edge_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        let right_pline_ = await exec_publics_Zone_node_9ylsi2n679i_extend_($p, right_edge_long_, 500);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(right_pline_, `name`, "right_pline");

        mfn.attrib.Set(right_pline_, `type`, right_type_);

        mfn.attrib.Set(right_pline_, `site_edge`, right_edge_long_);
      }
    }


    async function exec_publics_Zone_node_aj9qrn5t0cs_getLongestEdge_($p, edge_) {

      let road_ = mfn.attrib.Get(edge_, 'road');

      let edges_same_type_ = [edge_];

      let this_edge_ = edge_;

      let this_road_ = road_;

      while (this_road_ == road_) {

        this_edge_ = await exec_publics_Zone_node_aj9qrn5t0cs_nextEdge_($p, this_edge_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        this_road_ = mfn.attrib.Get(this_edge_, 'road');

        if (this_road_ == road_) {

          mfn.list.Add(edges_same_type_, this_edge_, 'to_end');
        }
      }

      this_edge_ = edge_;

      this_road_ = road_;

      while (this_road_ == road_) {

        this_edge_ = await exec_publics_Zone_node_aj9qrn5t0cs_prevEdge_($p, this_edge_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        this_road_ = mfn.attrib.Get(this_edge_, 'road');

        if (this_road_ == road_) {

          mfn.list.Add(edges_same_type_, this_edge_, 'to_end');
        }
      }

      edges_same_type_ = edges_same_type_;

      let lengths_ = mfn.calc.Length(edges_same_type_);

      return ifn.listSort(edges_same_type_, lengths_)[pythonList(-1, ifn.listSort(edges_same_type_, lengths_).length)];
    }


    async function exec_publics_Zone_node_aj9qrn5t0cs_getSiteLongestEdges_($p, site_, road_descr_) {

      let edges_ = [];

      for (let edge_ of mfn.query.Get('_e', site_)) {

        if (road_descr_ == mfn.attrib.Get(edge_, 'road')) {

          mfn.list.Add(edges_, edge_, 'to_end');
        }
      }

      if (ifn.len(edges_) == 0) {

        return [];
      }

      let posis_ = [mfn.query.Get('ps', edges_[pythonList(0, edges_.length)])];

      for (let edge_ of edges_.slice(1)) {

        let start_end_ = mfn.query.Get('ps', edge_);

        if (start_end_[pythonList(0, start_end_.length)] == posis_[pythonList(-1, posis_.length)][pythonList(-1, posis_[pythonList(-1, posis_.length)].length)]) {

          mfn.list.Add(posis_[pythonList(-1, posis_.length)], start_end_[pythonList(1, start_end_.length)], 'to_end');
        } else {
          if (start_end_[pythonList(1, start_end_.length)] == posis_[pythonList(-1, posis_.length)][pythonList(0, posis_[pythonList(-1, posis_.length)].length)]) {

            mfn.list.Add(posis_[pythonList(-1, posis_.length)], start_end_[pythonList(0, start_end_.length)], 'to_start');
          } else {

            mfn.list.Add(posis_, start_end_, 'to_end');
          }
        }
      }

      let site_plines_ = mfn.make.Polyline(posis_, 'open');

      return site_plines_;
    }


    async function exec_publics_Zone_node_aj9qrn5t0cs_extend_($p, pline_, dist_) {

      let edges_ = mfn.query.Get('_e', pline_);

      let lengths_ = mfn.calc.Length(edges_);

      let longest_ = ifn.listSort(edges_, lengths_)[pythonList(-1, ifn.listSort(edges_, lengths_).length)];

      let xyzs_ = mfn.attrib.Get(mfn.query.Get('ps', longest_), 'xyz');

      let vec_ = ifn.vecSetLen(ifn.vecFromTo(xyzs_[pythonList(0, xyzs_.length)], xyzs_[pythonList(1, xyzs_.length)]), dist_);

      let xyz0_ = ifn.vecSub(xyzs_[pythonList(0, xyzs_.length)], vec_);

      let xyz1_ = ifn.vecAdd(xyzs_[pythonList(1, xyzs_.length)], vec_);

      let posis_ = mfn.make.Position([xyz0_, xyz1_]);

      let long_pline_ = mfn.make.Polyline(posis_, 'open');

      return long_pline_;
    }


    async function exec_publics_Zone_node_aj9qrn5t0cs_offset_($p, sites_, roads_, dist_) {

      let off1_ = mfn.poly2d.OffsetMitre(roads_, dist_, dist_, 'square_end');

      let plines1_ = mfn.make.Polyline(off1_, 'close');

      let plines2_ = mfn.poly2d.Boolean(plines1_, sites_, 'intersect');

      mfn.edit.Delete([off1_, plines1_], 'delete_selected');

      return plines2_;
    }


    async function exec_publics_Zone_node_aj9qrn5t0cs_expandEdge_($p, pline_, dist_) {

      let pline_long_ = await exec_publics_Zone_node_aj9qrn5t0cs_extend_($p, pline_, 100);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let pgon_ = mfn.poly2d.OffsetMitre(pline_long_, dist_, dist_, 'square_end');

      return pgon_;
    }


    async function exec_publics_Zone_node_aj9qrn5t0cs_siteTrimPline_($p, site_, pline_) {

      let pline2_ = mfn.poly2d.Boolean(pline_, site_, 'intersect');

      mfn.edit.Delete(pline_, 'delete_selected');

      return pline2_;
    }


    async function exec_publics_Zone_node_aj9qrn5t0cs_getEdge_($p, wire_, posi_) {

      let posis_ = mfn.query.Get('ps', wire_);

      let edges_ = mfn.query.Get('_e', wire_);

      let idx_ = ifn.listFind(posis_, posi_);

      if (idx_ == -1) {

        return null;
      }

      return edges_[pythonList(idx_, edges_.length)];
    }


    async function exec_publics_Zone_node_aj9qrn5t0cs_nextEdge_($p, edge_) {

      let wire_ = mfn.query.Get('_w', edge_);

      let edges_ = mfn.query.Get('_e', wire_);

      let idx_ = ifn.listFind(edges_, edge_);

      let next_idx_ = idx_ + 1;

      if (next_idx_ == ifn.len(edges_)) {

        return edges_[pythonList(0, edges_.length)];
      }

      return edges_[pythonList(next_idx_, edges_.length)];
    }


    async function exec_publics_Zone_node_aj9qrn5t0cs_prevEdge_($p, edge_) {

      let wire_ = mfn.query.Get('_w', edge_);

      let edges_ = mfn.query.Get('_e', wire_);

      let idx_ = ifn.listFind(edges_, edge_);

      let prev_idx_ = idx_ - 1;

      if (prev_idx_ < 0) {

        return edges_[pythonList(-1, edges_.length)];
      }

      return edges_[pythonList(prev_idx_, edges_.length)];
    }


    async function exec_publics_Zone_node_aj9qrn5t0cs_sortAlongX_($p, plines_) {

      let cens_ = mfn.calc.Centroid(plines_, 'ps_average');

      let x_vals_ = [];

      for (let cen_ of cens_) {

        mfn.list.Add(x_vals_, cen_[pythonList(0, cen_.length)], 'to_end');
      }

      let sorted_ = ifn.listSort(plines_, x_vals_);

      let sorted_plines_ = [];

      for (let x_pline_ of sorted_) {

        mfn.list.Add(sorted_plines_, x_pline_[pythonList(1, x_pline_.length)], 'to_end');
      }

      return sorted_plines_;
    }


    async function exec_publics_Zone_node_aj9qrn5t0cs_closestEdge_($p, edges_, xyz_) {

      let min_d_ = Infinity;

      let closest_edge_ = null;

      if (xyz_[pythonList(0, xyz_.length)] > 1290) {
      }

      for (let edge_ of edges_) {

        let r_ = mfn.calc.Ray(edge_);

        let d_ = ifn.distance(xyz_, ifn.project(xyz_, r_, 0));

        if (d_ < min_d_) {

          min_d_ = d_;

          closest_edge_ = edge_;
        }

        if (d_ < 0.001) {

          break;
        }
      }

      return closest_edge_;
    }


    async function exec_publics_Zone_node_aj9qrn5t0cs_transferEdgeAttribs_($p, from_edges_, to_edges_) {

      for (let to_edge_ of to_edges_) {

        let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

        let from_edge_ = await exec_publics_Zone_node_aj9qrn5t0cs_closestEdge_($p, from_edges_, cen_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(to_edge_, `road`, mfn.attrib.Get(from_edge_, 'road'));
      }
    }

    async function exec_publics_Zone_node_aj9qrn5t0cs($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: zones', '__null__')
      }


      let site1_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', 'main_site')[pythonList(0, 'main_site'.length)];

      let art_pline_ = mfn.query.Filter(mfn.query.Get('pl', null), ['name', null], '==', "art_pline");

      let left_pline_ = mfn.query.Filter(mfn.query.Get('pl', null), ['name', null], '==', "left_pline");

      let right_pline_ = mfn.query.Filter(mfn.query.Get('pl', null), ['name', null], '==', "right_pline");

      let road_art_off_ = await exec_publics_Zone_node_aj9qrn5t0cs_offset_($p, site1_, art_pline_, BLK_ART_D_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      mfn.attrib.Set(road_art_off_, `name`, "road_loc_art_off");

      mfn.attrib.Set(road_art_off_, `type`, "road_loc_para");

      mfn.attrib.Set(mfn.query.Get('_e', road_art_off_), `road`, "road_loc_art_off");

      let left_off_ = [];

      let left_ = [];

      let road_left_off_ = [];

      if (ifn.len(left_pline_) > 0) {

        mfn.modify.Offset(left_pline_, -BLK_SEC_D_);

        let vec_ = mfn.calc.Vector(mfn.query.Get('_e', left_pline_)[pythonList(0, mfn.query.Get('_e', left_pline_).length)]);

        let ex_vec_ = ifn.vecSetLen([vec_[pythonList(1, vec_.length)], -vec_[pythonList(0, vec_.length)], 0], 300);

        left_off_ = mfn.make.Extrude(left_pline_, ex_vec_, 1, 'quads');

        mfn.edit.Reverse(left_off_);

        left_ = mfn.poly2d.Boolean(site1_, left_off_, 'intersect');

        road_left_off_ = await exec_publics_Zone_node_aj9qrn5t0cs_siteTrimPline_($p, site1_, left_pline_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(road_left_off_, `name`, "road_loc_left_off");

        mfn.attrib.Set(road_left_off_, `type`, "road_loc_perp");

        mfn.attrib.Set(mfn.query.Get('_e', road_left_off_), `road`, "road_loc_left_off");
      }

      let right_off_ = [];

      let right_ = [];

      let road_right_off_ = [];

      if (ifn.len(right_pline_) > 0) {

        mfn.modify.Offset(right_pline_, -BLK_SEC_D_);

        let vec_ = mfn.calc.Vector(mfn.query.Get('_e', right_pline_)[pythonList(0, mfn.query.Get('_e', right_pline_).length)]);

        let ex_vec_ = ifn.vecSetLen([vec_[pythonList(1, vec_.length)], -vec_[pythonList(0, vec_.length)], 0], 300);

        right_off_ = mfn.make.Extrude(right_pline_, ex_vec_, 1, 'quads');

        mfn.edit.Reverse(right_off_);

        right_ = mfn.poly2d.Boolean(site1_, right_off_, 'intersect');

        road_right_off_ = await exec_publics_Zone_node_aj9qrn5t0cs_siteTrimPline_($p, site1_, right_pline_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(road_right_off_, `name`, "road_loc_right_off");

        mfn.attrib.Set(road_right_off_, `type`, "road_loc_perp");

        mfn.attrib.Set(mfn.query.Get('_e', road_right_off_), `road`, "road_loc_right_off");
      }

      let dist_ = BLK_ART_D_;

      let art_off_ = mfn.poly2d.OffsetMitre(art_pline_, dist_, dist_, 'square_end');

      let off_grid_ = mfn.poly2d.Boolean(site1_, [left_off_, right_off_, art_off_], 'difference');

      let art_ = mfn.poly2d.Boolean(site1_, art_off_, 'intersect');

      let corners_ = mfn.poly2d.Boolean(art_, [left_off_, right_off_], 'intersect');

      let left2_ = mfn.poly2d.Boolean(left_, art_off_, 'difference');

      let right2_ = mfn.poly2d.Boolean(right_, art_off_, 'difference');

      let art2_ = mfn.poly2d.Boolean(art_, [left_off_, right_off_], 'difference');

      mfn.attrib.Set(off_grid_, `name`, "zone_mid");

      mfn.attrib.Set(art2_, `name`, "zone_art");

      mfn.attrib.Set(left2_, `name`, "zone_left");

      mfn.attrib.Set(right2_, `name`, "zone_right");

      mfn.attrib.Set(corners_, `name`, "zone_corner");

      let roads_ = [road_art_off_, road_left_off_, road_right_off_];

      let zones_ = ifn.listFlat([off_grid_, corners_, left2_, right2_, art2_]);

      let from_edges_ = ifn.listFlat([mfn.query.Get('_e', roads_), mfn.query.Get('_e', site1_)]);

      let to_edges_ = ifn.listFlat(mfn.query.Get('_e', zones_));

      await exec_publics_Zone_node_aj9qrn5t0cs_transferEdgeAttribs_($p, from_edges_, to_edges_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      mfn.edit.Delete([roads_, zones_], 'keep_selected');
    }


    async function exec_publics_Zone_node_z5vxdi5boqk_siteTrimPline_($p, site_, plines_) {

      let new_plines_ = [];

      for (let pline_ of plines_) {

        let pline2_ = mfn.poly2d.Boolean(pline_, site_, 'intersect');

        mfn.edit.Delete(pline_, 'delete_selected');

        mfn.list.Add(new_plines_, pline2_, 'to_end');
      }

      return new_plines_;
    }

    async function exec_publics_Zone_node_z5vxdi5boqk($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: art loc roads', '__null__')
      }


      let zone_ = mfn.query.Filter(mfn.query.Get('pg', null), ['name', null], '==', "zone_art");

      mfn.edit.Delete(zone_, 'keep_selected');

      let bb_ = mfn.calc.BBox(zone_);

      let num_x_ = ifn.round(bb_[pythonList(3, bb_.length)][pythonList(0, bb_[pythonList(3, bb_.length)].length)] / BLK_ART_W_) + 1;

      let x_plines_ = [];

      if (num_x_ > 2) {

        let posis_ = mfn.pattern.Grid(bb_[pythonList(0, bb_.length)], [bb_[pythonList(3, bb_.length)][pythonList(0, bb_[pythonList(3, bb_.length)].length)], bb_[pythonList(3, bb_.length)][pythonList(1, bb_[pythonList(3, bb_.length)].length)]], [num_x_, 2], 'columns');

        x_plines_ = mfn.make.Polyline(posis_, 'open');

        x_plines_ = x_plines_.slice(1, -1);

        x_plines_ = await exec_publics_Zone_node_z5vxdi5boqk_siteTrimPline_($p, zone_, x_plines_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(x_plines_, `type`, "road_loc_perp");
      }

      mfn.edit.Delete(x_plines_, 'keep_selected');
    }


    async function exec_publics_Zone_node_03lym4l18kqj_siteTrimPline_($p, site_, plines_) {

      let new_plines_ = [];

      for (let pline_ of plines_) {

        let pline2_ = mfn.poly2d.Boolean(pline_, site_, 'intersect');

        mfn.edit.Delete(pline_, 'delete_selected');

        mfn.list.Add(new_plines_, pline2_, 'to_end');
      }

      return new_plines_;
    }


    async function exec_publics_Zone_node_03lym4l18kqj_createRoads_($p, bb_all_, zone_, edges_off_) {

      let posis_left_off_ = mfn.query.Get('ps', edges_off_);

      let xyz0_ = mfn.attrib.Get(mfn.query.Get('ps', edges_off_[pythonList(0, edges_off_.length)])[pythonList(0, mfn.query.Get('ps', edges_off_[pythonList(0, edges_off_.length)]).length)], 'xyz');

      let xyz1_ = mfn.attrib.Get(mfn.query.Get('ps', edges_off_[pythonList(-1, edges_off_.length)])[pythonList(1, mfn.query.Get('ps', edges_off_[pythonList(-1, edges_off_.length)]).length)], 'xyz');

      let vec_ = ifn.vecFromTo(xyz0_, xyz1_);

      let bb_zone_ = mfn.calc.BBox(zone_);

      let dist_y_ = bb_zone_[pythonList(2, bb_zone_.length)][pythonList(1, bb_zone_[pythonList(2, bb_zone_.length)].length)] - bb_all_[pythonList(1, bb_all_.length)][pythonList(1, bb_all_[pythonList(1, bb_all_.length)].length)];

      let num_blocks_y_ = ifn.floor(dist_y_ / BLK_SEC_W_) + 1;

      let size_y_ = num_blocks_y_ * BLK_SEC_W_;

      let new_plines_ = [];

      if (num_blocks_y_ > 1) {

        for (let i_ of ifn.range(1, num_blocks_y_)) {

          let y_ = bb_all_[pythonList(1, bb_all_.length)][pythonList(1, bb_all_[pythonList(1, bb_all_.length)].length)] + (BLK_SEC_W_ * i_);

          let posi1_ = mfn.make.Position([bb_zone_[pythonList(1, bb_zone_.length)][pythonList(0, bb_zone_[pythonList(1, bb_zone_.length)].length)], y_, 0]);

          let posi2_ = mfn.make.Position([bb_zone_[pythonList(2, bb_zone_.length)][pythonList(0, bb_zone_[pythonList(2, bb_zone_.length)].length)], y_, 0]);

          let new_pline_ = mfn.make.Polyline([posi1_, posi2_], 'open');

          mfn.list.Add(new_plines_, new_pline_, 'to_end');
        }

        new_plines_ = await exec_publics_Zone_node_03lym4l18kqj_siteTrimPline_($p, zone_, new_plines_);
        if ($p.terminated) {
          return mfn.getModel();
        }
      }

      return new_plines_;
    }

    async function exec_publics_Zone_node_03lym4l18kqj($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: left and right loc roads', '__null__')
      }


      let zone_left_ = mfn.query.Filter(mfn.query.Get('pg', null), ['name', null], '==', "zone_left");

      let zone_mid_ = mfn.query.Filter(mfn.query.Get('pg', null), ['name', null], '==', "zone_mid");

      let zone_right_ = mfn.query.Filter(mfn.query.Get('pg', null), ['name', null], '==', "zone_right");

      let bb_all_ = mfn.calc.BBox(ifn.listFlat([zone_left_, zone_mid_, zone_right_]));

      mfn.edit.Delete([zone_left_, zone_right_], 'keep_selected');

      let edges_left_off_ = mfn.query.Filter(mfn.query.Get('_e', zone_left_), ['road', null], '==', "road_loc_left_off");

      let roads_left_ = [];

      if (ifn.len(edges_left_off_) > 0) {

        roads_left_ = await exec_publics_Zone_node_03lym4l18kqj_createRoads_($p, bb_all_, zone_left_, edges_left_off_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(roads_left_, `type`, "road_loc_para");
      }

      let edges_right_off_ = mfn.query.Filter(mfn.query.Get('_e', zone_right_), ['road', null], '==', "road_loc_right_off");

      let roads_right_ = [];

      if (ifn.len(edges_right_off_) > 0) {

        roads_right_ = await exec_publics_Zone_node_03lym4l18kqj_createRoads_($p, bb_all_, zone_right_, edges_right_off_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(roads_right_, `type`, "road_loc_para");
      }

      mfn.edit.Delete([roads_left_, roads_right_], 'keep_selected');
    }


    async function exec_publics_Zone_node_vu1s15ynip_siteTrimPline_($p, site_, plines_) {

      let new_plines_ = [];

      for (let pline_ of plines_) {

        let pline2_ = mfn.poly2d.Boolean(pline_, site_, 'intersect');

        mfn.edit.Delete(pline_, 'delete_selected');

        mfn.list.Add(new_plines_, pline2_, 'to_end');
      }

      return new_plines_;
    }

    async function exec_publics_Zone_node_vu1s15ynip($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: inner loc roads', '__null__')
      }


      let zone_left_ = mfn.query.Filter(mfn.query.Get('pg', null), ['name', null], '==', "zone_left");

      let zone_mid_ = mfn.query.Filter(mfn.query.Get('pg', null), ['name', null], '==', "zone_mid");

      let zone_right_ = mfn.query.Filter(mfn.query.Get('pg', null), ['name', null], '==', "zone_right");

      let bb_mid_ = mfn.calc.BBox(zone_mid_);

      let bb_all_ = mfn.calc.BBox(ifn.listFlat([zone_left_, zone_mid_, zone_right_]));

      mfn.edit.Delete(zone_mid_, 'keep_selected');

      let num_blks_x_ = ifn.round(bb_mid_[pythonList(3, bb_mid_.length)][pythonList(0, bb_mid_[pythonList(3, bb_mid_.length)].length)] / BLK_LOC_W_);

      let dist_y_ = bb_mid_[pythonList(2, bb_mid_.length)][pythonList(1, bb_mid_[pythonList(2, bb_mid_.length)].length)] - bb_all_[pythonList(1, bb_all_.length)][pythonList(1, bb_all_[pythonList(1, bb_all_.length)].length)];

      let num_blks_y_ = ifn.round(dist_y_ / BLK_LOC_D_);

      let cen_y_ = bb_all_[pythonList(1, bb_all_.length)][pythonList(1, bb_all_[pythonList(1, bb_all_.length)].length)] + (num_blks_y_ * BLK_LOC_D_ * 0.5);

      let plines_x_ = [];

      let plines_y_ = [];

      if (num_blks_x_ > 1 && num_blks_y_ > 1) {

        let origin_ = [bb_mid_[pythonList(0, bb_mid_.length)][pythonList(0, bb_mid_[pythonList(0, bb_mid_.length)].length)], cen_y_, 0];

        let size_ = [bb_mid_[pythonList(3, bb_mid_.length)][pythonList(0, bb_mid_[pythonList(3, bb_mid_.length)].length)], num_blks_y_ * BLK_LOC_D_];

        let num_posis_ = [num_blks_x_ + 1, num_blks_y_ + 1];

        let posis_ = mfn.pattern.Grid(origin_, size_, num_posis_, 'rows');

        plines_x_ = mfn.make.Polyline(posis_, 'open');

        plines_y_ = mfn.make.Polyline(ifn.listZip(posis_), 'open');

        plines_x_ = plines_x_.slice(1, -1);

        plines_x_ = await exec_publics_Zone_node_vu1s15ynip_siteTrimPline_($p, zone_mid_, plines_x_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(plines_x_, `type`, "road_loc_para_inn");

        plines_y_ = plines_y_.slice(1, -1);

        plines_y_ = await exec_publics_Zone_node_vu1s15ynip_siteTrimPline_($p, zone_mid_, plines_y_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(plines_y_, `type`, "road_loc_perp_inn");
      }

      mfn.edit.Delete([plines_x_, plines_y_], 'keep_selected');

      let var_ = mfn.edit.Fuse(mfn.query.Get('pl', null), 0.01);
    }


    async function exec_publics_Zone_node_ulnws3299id($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: local ring roads', '__null__')
      }


      mfn.edit.Delete(mfn.query.Get('pl', null), 'keep_selected');
    }


    async function exec_publics_Zone_node_93en2qi3yth($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: site', '__null__')
      }


      mfn.edit.Delete(mfn.query.Get('pg', null), 'keep_selected');
    }


    async function exec_publics_Zone_node_cpfgl91qjpc_touchingEdge_($p, edges_, xyz_) {

      for (let edge_ of edges_) {

        let r_ = mfn.calc.Ray(edge_);

        let d_ = ifn.distance(xyz_, ifn.project(xyz_, r_, 0));

        if (d_ < 0.01) {

          return edge_;

          break;
        }
      }

      return null;
    }


    async function exec_publics_Zone_node_cpfgl91qjpc_extend_($p, edge_, dist_) {

      let posis_ = mfn.query.Get('ps', edge_);

      let xyzs_ = mfn.attrib.Get(posis_, 'xyz');

      let vec_ = ifn.vecSetLen(ifn.vecFromTo(xyzs_[pythonList(0, xyzs_.length)], xyzs_[pythonList(1, xyzs_.length)]), dist_);

      if (dist_ > 0) {

        mfn.modify.Move(posis_[pythonList(1, posis_.length)], vec_);

        mfn.attrib.Set(posis_[pythonList(1, posis_.length)], `type`, "extended");
      } else {

        mfn.modify.Move(posis_[pythonList(0, posis_.length)], vec_);

        mfn.attrib.Set(posis_[pythonList(0, posis_.length)], `type`, "extended");
      }
    }


    async function exec_publics_Zone_node_cpfgl91qjpc_detachRoadEndsFromColdEdge_($p, site_edges_, roads_) {

      let roads_trimmed_ = [];

      for (let road_ of roads_) {

        let still_exists_ = await exec_publics_Zone_node_cpfgl91qjpc_detachFromColdEdge_($p, site_edges_, road_, true);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (still_exists_) {

          still_exists_ = await exec_publics_Zone_node_cpfgl91qjpc_detachFromColdEdge_($p, site_edges_, road_, false);
          if ($p.terminated) {
            return mfn.getModel();
          }

          if (still_exists_) {

            mfn.list.Add(roads_trimmed_, road_, 'to_end');
          }
        }
      }

      return roads_trimmed_;
    }


    async function exec_publics_Zone_node_cpfgl91qjpc_detachFromColdEdge_($p, site_edges_, road_, detach_start_) {

      let verts_ = mfn.query.Get('_v', road_);

      if (ifn.len(verts_) < 2) {

        return false;
      }

      let index_ = 0;

      if (!detach_start_) {

        index_ = -1;
      }

      let vert_ = verts_[pythonList(index_, verts_.length)];

      let posi_ = mfn.query.Get('ps', vert_)[pythonList(0, mfn.query.Get('ps', vert_).length)];

      let xyz_ = mfn.attrib.Get(posi_, 'xyz');

      let start_edge_ = await exec_publics_Zone_node_cpfgl91qjpc_touchingEdge_($p, site_edges_, xyz_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      if (start_edge_ != null) {

        let start_road_ = mfn.attrib.Get(start_edge_, 'road');

        if (start_road_ == "cold") {

          let var_ = mfn.edit.Weld([vert_], 'break_weld');

          mfn.edit.Delete(mfn.query.Get('ps', vert_), 'delete_selected');

          let edges_ = mfn.query.Get('_e', road_);

          if (ifn.len(edges_)) {

            if (detach_start_) {

              await exec_publics_Zone_node_cpfgl91qjpc_remDanglingEdges_($p, edges_[pythonList(0, edges_.length)]);
              if ($p.terminated) {
                return mfn.getModel();
              }
            } else {

              await exec_publics_Zone_node_cpfgl91qjpc_remDanglingEdges_($p, edges_[pythonList(-1, edges_.length)]);
              if ($p.terminated) {
                return mfn.getModel();
              }
            }
          }
        }
      }

      let exists_ = mfn.query.Type(road_, 'exists');

      return exists_;
    }


    async function exec_publics_Zone_node_cpfgl91qjpc_remDanglingEdges_($p, edge_) {

      let road_ = mfn.query.Get('pl', edge_);

      let edges_ = mfn.query.Get('_e', road_);

      let index_ = 0;

      if (edges_[pythonList(-1, edges_.length)] == edge_) {

        mfn.list.Sort(edges_, 'reverse');

        index_ = 1;
      }

      for (edge_ of edges_) {

        let edge_posi_ = mfn.query.Get('ps', edge_)[pythonList(index_, mfn.query.Get('ps', edge_).length)];

        if (ifn.len(mfn.query.Get('_e', edge_posi_)) == 1) {

          mfn.edit.Delete(edge_posi_, 'delete_selected');
        } else {

          break;
        }
      }
    }


    async function exec_publics_Zone_node_cpfgl91qjpc_extendBeyondBoundaryRoad_($p, site_edges_, roads_) {

      for (let road_ of roads_) {

        let posis_ = mfn.query.Get('ps', road_);

        let edges_ = mfn.query.Get('_e', road_);

        let start_edge_ = await exec_publics_Zone_node_cpfgl91qjpc_touchingEdge_($p, site_edges_, mfn.attrib.Get(posis_[pythonList(0, posis_.length)], 'xyz'));
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (start_edge_ != null) {

          let start_road_ = mfn.attrib.Get(start_edge_, 'road');

          if (ifn.listFind(["road_art", "road_sec", "road_ter"], start_road_) != -1) {

            await exec_publics_Zone_node_cpfgl91qjpc_extend_($p, edges_[pythonList(0, edges_.length)], -20);
            if ($p.terminated) {
              return mfn.getModel();
            }
          }
        }

        let end_edge_ = await exec_publics_Zone_node_cpfgl91qjpc_touchingEdge_($p, site_edges_, mfn.attrib.Get(posis_[pythonList(-1, posis_.length)], 'xyz'));
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (end_edge_ != null) {

          let end_road_ = mfn.attrib.Get(end_edge_, 'road');

          if (ifn.listFind(["road_art", "road_sec", "road_ter"], end_road_) != -1) {

            await exec_publics_Zone_node_cpfgl91qjpc_extend_($p, edges_[pythonList(-1, edges_.length)], 20);
            if ($p.terminated) {
              return mfn.getModel();
            }
          }
        }
      }
    }


    async function exec_publics_Zone_node_cpfgl91qjpc_angDot_($p, vert_) {

      let edges_ = mfn.query.Get('_e', vert_);

      let vecs_ = mfn.calc.Vector(edges_);

      let vec0_ = ifn.vecNorm(vecs_[pythonList(0, vecs_.length)]);

      let vec1_ = ifn.vecNorm(vecs_[pythonList(1, vecs_.length)]);

      return ifn.vecDot(vec0_, vec1_);
    }

    async function exec_publics_Zone_node_cpfgl91qjpc($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: roads', '__null__')
      }


      let site1_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', 'main_site')[pythonList(0, 'main_site'.length)];

      let roads_ = mfn.query.Get('pl', null);

      let new_roads_ = mfn.poly2d.Stitch(roads_, 0.1);

      mfn.edit.Delete(roads_, 'delete_selected');

      let site_edges_ = mfn.query.Get('_e', site1_);

      new_roads_ = await exec_publics_Zone_node_cpfgl91qjpc_detachRoadEndsFromColdEdge_($p, site_edges_, new_roads_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      new_roads_ = mfn.query.Get('pl', new_roads_);

      await exec_publics_Zone_node_cpfgl91qjpc_extendBeyondBoundaryRoad_($p, site_edges_, new_roads_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      new_roads_ = mfn.query.Get('pl', new_roads_);

      let var_ = mfn.edit.Fuse(new_roads_, 5);

      for (let posi_ of mfn.query.Get('ps', new_roads_)) {

        if (ifn.len(mfn.query.Get('_e', posi_)) == 1 && mfn.attrib.Get(posi_, 'type') != "extended") {

          await exec_publics_Zone_node_cpfgl91qjpc_remDanglingEdges_($p, mfn.query.Get('_e', posi_));
          if ($p.terminated) {
            return mfn.getModel();
          }
        }
      }

      new_roads_ = mfn.query.Get('pl', new_roads_);

      for (let posi_ of mfn.query.Get('ps', new_roads_)) {

        let verts_ = mfn.query.Get('_v', posi_);

        let edges_ = mfn.query.Get('_e', posi_);

        if (ifn.len(edges_) == 2 && ifn.len(verts_) == 1) {

          let dot_ = await exec_publics_Zone_node_cpfgl91qjpc_angDot_($p, verts_[pythonList(0, verts_.length)]);
          if ($p.terminated) {
            return mfn.getModel();
          }

          if (ifn.abs(dot_) > 0.99) {

            mfn.edit.Delete(posi_, 'delete_selected');
          }
        }
      }

      new_roads_ = mfn.query.Get('pl', new_roads_);

      mfn.edit.Delete([site1_, new_roads_], 'keep_selected');
    }


    async function exec_publics_Zone_node_8joyldrk8dk_transferEdgeAttribs_($p, from_edges_, to_edges_) {

      for (let to_edge_ of to_edges_) {

        let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

        let from_edge_ = await exec_publics_Zone_node_8joyldrk8dk_touchingEdge_($p, from_edges_, cen_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (from_edge_ != null) {

          mfn.attrib.Set(to_edge_, `road`, mfn.attrib.Get(from_edge_, 'road'));
        } else {

          mfn.attrib.Set(to_edge_, `road`, "road_loc");
        }
      }
    }


    async function exec_publics_Zone_node_8joyldrk8dk_touchingEdge_($p, edges_, xyz_) {

      for (let edge_ of edges_) {

        let r_ = mfn.calc.Ray(edge_);

        let d_ = ifn.distance(xyz_, ifn.project(xyz_, r_, 0));

        if (d_ < 1) {

          return edge_;
        }
      }
    }


    async function exec_publics_Zone_node_8joyldrk8dk_joinPlines_($p, pline0_, pline1_) {

      let posis0_ = mfn.query.Get('ps', pline0_);

      let posis1_ = mfn.query.Get('ps', pline1_);

      let posis_ = null;

      if (posis0_[pythonList(0, posis0_.length)] == posis1_[pythonList(0, posis1_.length)]) {

        posis_ = ifn.listJoin(ifn.listRev(posis0_), posis1_.slice(1));
      } else {
        if (posis0_[pythonList(-1, posis0_.length)] == posis1_[pythonList(-1, posis1_.length)]) {

          posis_ = ifn.listJoin(posis0_, ifn.listRev(posis1_).slice(1));
        } else {
          if (posis0_[pythonList(-1, posis0_.length)] == posis1_[pythonList(0, posis1_.length)]) {

            posis_ = ifn.listJoin(posis0_, posis1_.slice(1));
          } else {
            if (posis0_[pythonList(0, posis0_.length)] == posis1_[pythonList(-1, posis1_.length)]) {

              posis_ = ifn.listJoin(ifn.listRev(posis0_), ifn.listRev(posis1_).slice(1));
            }
          }
        }
      }

      let cc_ = posis_;

      let new_pline_ = mfn.make.Polyline(posis_, 'open');

      mfn.edit.Delete([pline0_, pline1_], 'delete_selected');

      return new_pline_;
    }


    async function exec_publics_Zone_node_8joyldrk8dk__cleanPgonsEdge_($p, pgons_) {

      for (let pgon_ of pgons_) {

        let del_posis_ = [];

        for (let edge_ of mfn.query.Get('_e', pgon_)) {

          let length_ = mfn.calc.Length(edge_);

          if (length_ < 1) {

            let posis_ = mfn.query.Get('ps', edge_);

            mfn.list.Add(del_posis_, posis_[pythonList(0, posis_.length)], 'to_end');
          }
        }

        mfn.edit.Delete(del_posis_, 'delete_selected');
      }

      return mfn.query.Get('pg', pgons_);
    }


    async function exec_publics_Zone_node_8joyldrk8dk__cleanPgonsAng_($p, pgons_) {

      for (let pgon_ of pgons_) {

        let del_posis_ = [];

        for (let vert_ of mfn.query.Get('_v', pgon_)) {

          let dot_ = await exec_publics_Zone_node_8joyldrk8dk__angDot_($p, vert_);
          if ($p.terminated) {
            return mfn.getModel();
          }

          if (ifn.abs(dot_) > 0.9999) {

            mfn.list.Add(del_posis_, mfn.query.Get('ps', vert_), 'to_end');
          }
        }

        mfn.edit.Delete(del_posis_, 'delete_selected');
      }

      return mfn.query.Get('pg', pgons_);
    }


    async function exec_publics_Zone_node_8joyldrk8dk__angDot_($p, vert_) {

      let edges_ = mfn.query.Get('_e', vert_);

      let vecs_ = mfn.calc.Vector(edges_);

      let vec0_ = ifn.vecNorm(vecs_[pythonList(0, vecs_.length)]);

      let vec1_ = ifn.vecNorm(vecs_[pythonList(1, vecs_.length)]);

      return ifn.vecDot(vec0_, vec1_);
    }


    async function exec_publics_Zone_node_8joyldrk8dk_cleanPgons_($p, pgons_) {

      pgons_ = await exec_publics_Zone_node_8joyldrk8dk__cleanPgonsEdge_($p, pgons_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      pgons_ = await exec_publics_Zone_node_8joyldrk8dk__cleanPgonsAng_($p, pgons_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      return pgons_;
    }

    async function exec_publics_Zone_node_8joyldrk8dk($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: blocks', '__null__')
      }


      let site1_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', 'main_site')[pythonList(0, 'main_site'.length)];

      let roads_ = mfn.query.Get('pl', null);

      let road_pgons_ = [];

      let blocks_ = [];

      if (ifn.len(roads_)) {

        let bool_roads_ = mfn.make.Copy(roads_, null);

        for (let posi_ of mfn.query.Get('ps', bool_roads_)) {

          let edges_ = mfn.query.Get('_e', posi_);

          let plines_ = mfn.query.Get('pl', posi_);

          if (ifn.len(plines_) == 2 && ifn.len(edges_) == 2) {

            let new_pline_ = await exec_publics_Zone_node_8joyldrk8dk_joinPlines_($p, plines_[pythonList(0, plines_.length)], plines_[pythonList(1, plines_.length)]);
            if ($p.terminated) {
              return mfn.getModel();
            }

            mfn.list.Add(bool_roads_, new_pline_, 'to_end');
          }
        }

        bool_roads_ = mfn.query.Get('pl', bool_roads_);

        road_pgons_ = mfn.poly2d.OffsetMitre(mfn.query.Get('pl', bool_roads_), ROAD_LOC_W_ / 2, 100, 'butt_end');

        road_pgons_ = await exec_publics_Zone_node_8joyldrk8dk_cleanPgons_($p, road_pgons_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        blocks_ = mfn.poly2d.Boolean(site1_, road_pgons_, 'difference');

        blocks_ = await exec_publics_Zone_node_8joyldrk8dk_cleanPgons_($p, blocks_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        road_pgons_ = mfn.poly2d.Boolean(road_pgons_, site1_, 'intersect');

        road_pgons_ = mfn.poly2d.Union(road_pgons_);

        let site_edges_ = mfn.query.Get('_e', site1_);

        for (let block_ of blocks_) {

          await exec_publics_Zone_node_8joyldrk8dk_transferEdgeAttribs_($p, site_edges_, mfn.query.Get('_e', block_));
          if ($p.terminated) {
            return mfn.getModel();
          }
        }
      } else {

        blocks_ = mfn.make.Copy(site1_, null);

        let site_edges_ = mfn.query.Get('_e', site1_);

        await exec_publics_Zone_node_8joyldrk8dk_transferEdgeAttribs_($p, site_edges_, mfn.query.Get('_e', blocks_));
        if ($p.terminated) {
          return mfn.getModel();
        }
      }

      mfn.edit.Delete([roads_, blocks_, road_pgons_], 'keep_selected');

      mfn.attrib.Set(blocks_, `type`, "block");

      mfn.attrib.Set(road_pgons_, `type`, "road_loc");
    }


    async function exec_publics_Zone_node_8tqxj3gvnfc($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: undo rot', '__null__')
      }


      let ents_ = [mfn.query.Get('pl', null), mfn.query.Get('pg', null)];

      mfn.modify.XForm(ents_, JSON.parse(JSON.stringify(ifn.XY)), mfn.attrib.Get(null, 'pln'));
    }


    async function exec_publics_Zone_node_9bhi1flny5t($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: End', '__null__')
      }


      return null;
    }

    var merged;
    let ssid_exec_publics_Zone_node_8ylatkqudwy = mfn.model.snapshotGetActive();

    let ssid_exec_publics_Zone_node_xlp2w5q5fr = mfn.model.snapshotNext([ssid_exec_publics_Zone_node_8ylatkqudwy]);

    await exec_publics_Zone_node_xlp2w5q5fr($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_publics_Zone_node_9ylsi2n679i = mfn.model.snapshotNext([ssid_exec_publics_Zone_node_xlp2w5q5fr]);

    await exec_publics_Zone_node_9ylsi2n679i($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_publics_Zone_node_aj9qrn5t0cs = ssid_exec_publics_Zone_node_9ylsi2n679i;

    await exec_publics_Zone_node_aj9qrn5t0cs($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_publics_Zone_node_z5vxdi5boqk = mfn.model.snapshotNext([ssid_exec_publics_Zone_node_aj9qrn5t0cs]);

    await exec_publics_Zone_node_z5vxdi5boqk($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_publics_Zone_node_03lym4l18kqj = mfn.model.snapshotNext([ssid_exec_publics_Zone_node_aj9qrn5t0cs]);

    await exec_publics_Zone_node_03lym4l18kqj($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_publics_Zone_node_vu1s15ynip = mfn.model.snapshotNext([ssid_exec_publics_Zone_node_aj9qrn5t0cs]);

    await exec_publics_Zone_node_vu1s15ynip($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_publics_Zone_node_ulnws3299id = mfn.model.snapshotNext([ssid_exec_publics_Zone_node_aj9qrn5t0cs]);

    await exec_publics_Zone_node_ulnws3299id($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_publics_Zone_node_93en2qi3yth = mfn.model.snapshotNext([ssid_exec_publics_Zone_node_xlp2w5q5fr]);

    await exec_publics_Zone_node_93en2qi3yth($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_publics_Zone_node_cpfgl91qjpc = mfn.model.snapshotNext([ssid_exec_publics_Zone_node_z5vxdi5boqk, ssid_exec_publics_Zone_node_03lym4l18kqj, ssid_exec_publics_Zone_node_vu1s15ynip, ssid_exec_publics_Zone_node_ulnws3299id, ssid_exec_publics_Zone_node_93en2qi3yth]);

    await exec_publics_Zone_node_cpfgl91qjpc($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_publics_Zone_node_8joyldrk8dk = ssid_exec_publics_Zone_node_cpfgl91qjpc;

    await exec_publics_Zone_node_8joyldrk8dk($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_publics_Zone_node_8tqxj3gvnfc = ssid_exec_publics_Zone_node_8joyldrk8dk;

    await exec_publics_Zone_node_8tqxj3gvnfc($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_publics_Zone_node_9bhi1flny5t = mfn.model.snapshotNext([ssid_exec_publics_Zone_node_8ylatkqudwy, ssid_exec_publics_Zone_node_8tqxj3gvnfc]);

    return await exec_publics_Zone_node_9bhi1flny5t($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);
  }

  async function exec_publics($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {

    async function exec_publics_node_huiuxvcjtij($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: import', '__null__')
      }


      let coll_city_ = await mfn.io.Import(IN_MODEL_, 'sim');
    }


    async function exec_publics_node_u4td3ir3jf_getCenBlk_($p, parts_, block_id_) {

      let ogblocks_ = mfn.query.Filter(mfn.query.Filter(parts_, ['block_type', null], '!=', 'art'), ['block_type', null], '!=', 'sec');

      let parts_cen_ = mfn.query.Filter(ogblocks_, ['block_id', null], '==', block_id_);

      let parts_other_ = mfn.query.Filter(ogblocks_, ['block_id', null], '!=', block_id_);

      return [parts_cen_, parts_other_];
    }


    async function exec_publics_node_u4td3ir3jf_processCenBlk_($p, parts_cen_, site_cen_, open_req_area_, amen_req_area_) {

      let total_req_area_ = open_req_area_ + amen_req_area_;

      parts_cen_ = parts_cen_;

      let reordered_parts_ = await exec_publics_node_u4td3ir3jf__reorderParts1_($p, parts_cen_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let result_ = await exec_publics_node_u4td3ir3jf__cenBlkOpenAndAmen_B_($p, reordered_parts_, site_cen_, open_req_area_, amen_req_area_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      return result_;
    }


    async function exec_publics_node_u4td3ir3jf_insertIntoBlks_($p, og_parts_edges_, loc_parts_, curr_area_, req_area_, avoid_nei_) {

      let result_ = await exec_publics_node_u4td3ir3jf__addOffGridPartsUntReq_($p, og_parts_edges_, curr_area_, req_area_, avoid_nei_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let added_parts_ = result_[pythonList(0, result_.length)];

      let overlap_parts_ = result_[pythonList(1, result_.length)];

      curr_area_ = result_[pythonList(2, result_.length)];

      let all_new_loc_parts_ = [];

      for (let loc_part_ of loc_parts_) {

        let block_id_ = mfn.attrib.Get(loc_part_, 'block_id');

        let olaps_ = mfn.query.Filter(overlap_parts_, ['block_id', null], '==', block_id_);

        let olaps_o_ = mfn.query.Filter(overlap_parts_, ['block_id', null], '==', block_id_ + "o");

        let olaps_r_ = mfn.query.Filter(overlap_parts_, ['block_id', null], '==', block_id_ + "r");

        olaps_ = ifn.listFlat([olaps_, olaps_o_, olaps_r_]);

        let isects_ = mfn.poly2d.Boolean(olaps_, loc_part_, 'intersect');

        mfn.list.Add(all_new_loc_parts_, isects_, 'extend_end');

        mfn.list.Add(added_parts_, isects_, 'extend_end');

        if (ifn.len(isects_) > 0) {

          let new_loc_parts_ = mfn.poly2d.Boolean(loc_part_, olaps_, 'difference');

          for (let new_loc_part_ of new_loc_parts_) {

            let area_ = mfn.calc.Area(new_loc_part_);

            if (area_ > 0.1) {

              mfn.list.Add(all_new_loc_parts_, new_loc_part_, 'extend_end');
            } else {

              mfn.edit.Delete(new_loc_part_, 'delete_selected');
            }
          }

          new_loc_parts_ = mfn.query.Get('pg', new_loc_parts_);

          await exec_publics_node_u4td3ir3jf__coptAttribs_($p, loc_part_, mfn.query.Get('pg', new_loc_parts_));
          if ($p.terminated) {
            return mfn.getModel();
          }

          await exec_publics_node_u4td3ir3jf__transferEdgeAttribs_($p, mfn.query.Get('_e', loc_part_), mfn.query.Get('_e', new_loc_parts_));
          if ($p.terminated) {
            return mfn.getModel();
          }

          await exec_publics_node_u4td3ir3jf__coptAttribs_($p, loc_part_, isects_);
          if ($p.terminated) {
            return mfn.getModel();
          }

          await exec_publics_node_u4td3ir3jf__transferEdgeAttribs_($p, mfn.query.Get('_e', loc_part_), mfn.query.Get('_e', isects_));
          if ($p.terminated) {
            return mfn.getModel();
          }

          await exec_publics_node_u4td3ir3jf__transferEdgeAttribsBtwTouchingParts_($p, ifn.listJoin(isects_, new_loc_parts_));
          if ($p.terminated) {
            return mfn.getModel();
          }

          mfn.edit.Delete(loc_part_, 'delete_selected');
        }
      }

      mfn.edit.Delete(overlap_parts_, 'delete_selected');

      return [added_parts_, all_new_loc_parts_, curr_area_];
    }


    async function exec_publics_node_u4td3ir3jf__addOffGridPartsUntReq_($p, og_parts_edges_, curr_area_, req_area_, avoid_nei_) {

      let added_parts_ = [];

      let overlap_parts_ = [];

      let nei_parts_ = [];

      if (curr_area_ > req_area_) {

        return [added_parts_, overlap_parts_, curr_area_];
      }

      for (let part_edge_ of og_parts_edges_) {

        let part_ = part_edge_[pythonList(0, part_edge_.length)];

        let type_ = mfn.attrib.Get(part_, 'type');

        let edge_ = part_edge_[pythonList(1, part_edge_.length)];

        let edge_type_ = mfn.attrib.Get(edge_, 'road');

        if (ifn.listHas(nei_parts_, part_)) {

          continue;
        }

        mfn.list.Add(added_parts_, part_, 'extend_end');

        let area_ = mfn.calc.Area(part_);

        curr_area_ = curr_area_ + area_;

        if (type_ == "concave_corner") {

          if (edge_) {

            let block_pgs_ = mfn.query.Filter(mfn.query.Get('pg', null), ['block_id', null], '==', mfn.attrib.Get(part_, 'block_id'));
            printFunc($p.console, `block_pgs`, block_pgs_);

            for (let pg_ of block_pgs_) {

              if (mfn.attrib.Get(pg_, 'type') != 'concave_corner') {

                let loc_overlap_ = mfn.make.Copy(pg_, null);

                mfn.list.Add(overlap_parts_, loc_overlap_, 'to_end');

                area_ = mfn.calc.Area(loc_overlap_);

                curr_area_ = curr_area_ + area_;
              }
            }
          }
        }

        if (type_ == "off_grid0" || type_ == "off_grid1") {

          if (edge_) {

            let edge_vec_ = mfn.calc.Vector(edge_);

            let olap_dist_ = PART_LOC_D_ + 2;

            if (edge_type_ == "sec") {

              olap_dist_ = PART_SEC_D_ + 2;
            }

            let perp_vec_ = ifn.vecSetLen([edge_vec_[pythonList(1, edge_vec_.length)], -edge_vec_[pythonList(0, edge_vec_.length)], 0], olap_dist_);

            let loc_overlap_ = mfn.make.Extrude(edge_, perp_vec_, 1, 'quads');

            let var_ = mfn.edit.Weld(loc_overlap_, 'break_weld');

            mfn.modify.Move(loc_overlap_, ifn.vecSetLen(perp_vec_, -1));

            mfn.attrib.Set(loc_overlap_, `block_id`, mfn.attrib.Get(part_, 'block_id'));

            mfn.list.Add(overlap_parts_, loc_overlap_, 'to_end');

            area_ = mfn.calc.Area(loc_overlap_);

            curr_area_ = curr_area_ + area_;
          }
        }

        if (type_ == "off_grid1") {

          let cluster_id_ = mfn.attrib.Get(part_, 'cluster_id');

          if (!ifn.isUndef(cluster_id_)) {

            let cluster_parts_ = mfn.query.Filter(mfn.query.Get('pg', null), ['cluster_id', null], '==', cluster_id_);

            mfn.list.Remove(cluster_parts_, part_, 'first_value');

            for (let back_ of cluster_parts_) {

              if (ifn.listHas(nei_parts_, back_)) {

                continue;
              }

              mfn.list.Add(added_parts_, back_, 'extend_end');

              area_ = mfn.calc.Area(part_);

              curr_area_ = curr_area_ + area_;
            }
          }
        }

        if (curr_area_ > req_area_) {

          break;
        }

        if (avoid_nei_) {

          for (let added_part_ of added_parts_) {

            let nei_ = mfn.query.Neighbor('pg', part_);

            nei_ = mfn.query.Filter(nei_, ['type', null], '==', mfn.attrib.Get(part_, 'type'));

            mfn.list.Add(nei_parts_, nei_, 'extend_end');
          }
        }
      }

      return [added_parts_, overlap_parts_, curr_area_];
    }


    async function exec_publics_node_u4td3ir3jf__getBackParts_($p, part_edge_, type_) {

      let edges_ = mfn.query.Get('_e', part_edge_[pythonList(0, part_edge_.length)]);

      let len_edges_ = ifn.len(edges_);

      let idx_ = ifn.listFind(edges_, part_edge_[pythonList(1, part_edge_.length)]);

      let edges2_ = [];

      for (let i_ of ifn.range(len_edges_)) {

        mfn.list.Add(edges2_, edges_[pythonList((idx_ + i_) % len_edges_, edges_.length)], 'to_end');
      }

      let back_ = edges2_.slice(2, -1);

      let back_edges_ = mfn.query.Edge(back_, 'touching');

      back_edges_ = ifn.listFlat(back_edges_);

      let back_pgons_ = mfn.query.Filter(mfn.query.Get('pg', back_edges_), ['type', null], '==', type_);

      back_pgons_ = ifn.listFlat(back_pgons_);

      return back_pgons_;
    }


    async function exec_publics_node_u4td3ir3jf__cenBlkAllOpen_($p, parts_, open_req_area_) {

      let areas_ = mfn.calc.Area(parts_);

      let curr_open_area_ = ifn.sum(areas_);

      await exec_publics_node_u4td3ir3jf__makeOpen_($p, parts_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let result_open_ = await exec_publics_node_u4td3ir3jf__printStr_($p, curr_open_area_, open_req_area_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      return curr_open_area_;
    }


    async function exec_publics_node_u4td3ir3jf__cenBlkOpenAndAmen_A_($p, parts_, site_cen_, open_req_area_, amen_req_area_) {

      let parts_cen_ = await exec_publics_node_u4td3ir3jf__sortPartsStrightDistToCen_($p, parts_, site_cen_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      parts_cen_ = ifn.listRev(parts_cen_);

      let reordered_parts_ = await exec_publics_node_u4td3ir3jf__reorderParts_($p, parts_cen_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let curr_open_area_ = 0;

      let curr_open_ = [];

      curr_open_area_ = await exec_publics_node_u4td3ir3jf__addPartsUntilReqArea_($p, curr_open_, reordered_parts_, curr_open_area_, open_req_area_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let neighbours_ = mfn.query.Neighbor('pg', curr_open_);

      let corners_ = await exec_publics_node_u4td3ir3jf__filterKeepCorners_($p, neighbours_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let curr_amen_area_ = await exec_publics_node_u4td3ir3jf__addParts_($p, curr_open_, corners_, curr_open_area_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_publics_node_u4td3ir3jf__makeOpen_($p, curr_open_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let parts_remainder_ = await exec_publics_node_u4td3ir3jf__filterRem_($p, reordered_parts_, curr_open_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_publics_node_u4td3ir3jf__makeAmen_($p, parts_remainder_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let areas_ = mfn.calc.Area(parts_remainder_);

      curr_amen_area_ = ifn.sum(areas_);

      let result_open_ = await exec_publics_node_u4td3ir3jf__printStr_($p, curr_open_area_, open_req_area_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let result_amen_ = await exec_publics_node_u4td3ir3jf__printStr_($p, curr_amen_area_, amen_req_area_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      return [curr_open_area_, curr_amen_area_];
    }


    async function exec_publics_node_u4td3ir3jf__cenBlkOpenAndAmen_B_($p, parts_, site_cen_, open_req_area_, amen_req_area_) {

      let part_open_ = [];

      let parts_leftover_ = [];

      for (let i_ of ifn.range(ifn.len(parts_))) {

        if ((i_ > 1) || (i_ > ifn.len(parts_) - 1)) {

          mfn.list.Add(parts_leftover_, parts_[pythonList(i_, parts_.length)][pythonList(0, parts_[pythonList(i_, parts_.length)].length)], 'to_end');

          for (let part_other_ of parts_[pythonList(i_, parts_.length)][pythonList(3, parts_[pythonList(i_, parts_.length)].length)]) {

            mfn.list.Add(parts_leftover_, part_other_, 'to_end');
          }
        } else {

          mfn.list.Add(part_open_, parts_[pythonList(i_, parts_.length)][pythonList(0, parts_[pythonList(i_, parts_.length)].length)], 'to_end');

          for (let part_other_ of parts_[pythonList(i_, parts_.length)][pythonList(3, parts_[pythonList(i_, parts_.length)].length)]) {

            mfn.list.Add(part_open_, part_other_, 'to_end');
          }
        }
      }

      let curr_open_area_ = 0;

      let curr_open_ = [];

      curr_open_area_ = await exec_publics_node_u4td3ir3jf__addPartsUntilReqArea_($p, curr_open_, part_open_, curr_open_area_, open_req_area_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_publics_node_u4td3ir3jf__makeOpen_($p, curr_open_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let curr_amen_area_ = 0;

      let curr_amen_ = [];

      let parts_remainder_ = await exec_publics_node_u4td3ir3jf__filterRem_($p, part_open_, curr_open_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let parts_for_amen_ = await exec_publics_node_u4td3ir3jf__filterThroughpart_($p, [parts_remainder_, parts_leftover_]);
      if ($p.terminated) {
        return mfn.getModel();
      }

      curr_amen_area_ = await exec_publics_node_u4td3ir3jf__addPartsUntilReqArea_($p, curr_amen_, parts_for_amen_, curr_amen_area_, amen_req_area_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_publics_node_u4td3ir3jf__makeAmen_($p, curr_amen_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let main_park_hull_ = mfn.poly2d.ConvexHull(curr_open_);

      return [curr_open_area_, curr_amen_area_, main_park_hull_];
    }


    async function exec_publics_node_u4td3ir3jf__cenBlkOpenAndAmen_($p, parts_, site_cen_, open_req_area_, amen_req_area_) {

      let parts_cen_ = await exec_publics_node_u4td3ir3jf__sortPartsStrightDistToCen_($p, parts_, site_cen_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      parts_cen_ = ifn.listRev(parts_cen_);

      let result_ = await exec_publics_node_u4td3ir3jf__divParts_($p, parts_cen_);
      if ($p.terminated) {
        return mfn.getModel();
      }
      printFunc($p.console, `result`, result_);

      let parts_cor_ = result_[pythonList(0, result_.length)];

      let parts_side_ = result_[pythonList(1, result_.length)].slice(0, -2);

      let curr_amen_area_ = 0;

      let curr_amen_ = [];

      curr_amen_area_ = await exec_publics_node_u4td3ir3jf__addPartsUntilReqArea_($p, curr_amen_, parts_side_, curr_amen_area_, amen_req_area_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let neighbours_ = mfn.query.Neighbor('pg', curr_amen_);

      let corners_ = await exec_publics_node_u4td3ir3jf__filterKeepCorners_($p, neighbours_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      curr_amen_area_ = await exec_publics_node_u4td3ir3jf__addParts_($p, curr_amen_, corners_, curr_amen_area_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_publics_node_u4td3ir3jf__makeAmen_($p, curr_amen_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let parts_remainder_ = await exec_publics_node_u4td3ir3jf__filterRem_($p, parts_, curr_amen_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_publics_node_u4td3ir3jf__makeOpen_($p, parts_remainder_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let areas_ = mfn.calc.Area(parts_remainder_);

      let curr_open_area_ = ifn.sum(areas_);

      let result_open_ = await exec_publics_node_u4td3ir3jf__printStr_($p, curr_open_area_, open_req_area_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let result_amen_ = await exec_publics_node_u4td3ir3jf__printStr_($p, curr_amen_area_, amen_req_area_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      return [curr_open_area_, curr_amen_area_];
    }


    async function exec_publics_node_u4td3ir3jf__cenBlkMixed_($p, parts_, site_cen_, open_req_area_, amen_req_area_, resi_req_area_) {

      let parts_cen_ = await exec_publics_node_u4td3ir3jf__sortPartsStrightDistToCen_($p, parts_, site_cen_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      parts_cen_ = ifn.listRev(parts_cen_);

      let result_ = await exec_publics_node_u4td3ir3jf__divParts_($p, parts_cen_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let parts_cor_ = result_[pythonList(0, result_.length)];

      let parts_side_ = result_[pythonList(1, result_.length)];

      let parts_side0_ = parts_side_[pythonList(-1, parts_side_.length)];

      let parts_side_remainder_ = parts_side_.slice(0, -1);

      let parts_og_ = result_[pythonList(2, result_.length)];

      let curr_amen_ = [];

      let curr_amen_area_ = 0;

      curr_amen_area_ = await exec_publics_node_u4td3ir3jf__addPartsUntilReqArea_($p, curr_amen_, parts_side_remainder_, curr_amen_area_, amen_req_area_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      parts_side_remainder_ = await exec_publics_node_u4td3ir3jf__filterRem_($p, parts_side_remainder_, curr_amen_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let neighbours_ = mfn.query.Neighbor('pg', curr_amen_);

      let amen_corners_ = await exec_publics_node_u4td3ir3jf__filterKeepCorners_($p, neighbours_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      curr_amen_area_ = await exec_publics_node_u4td3ir3jf__addParts_($p, curr_amen_, amen_corners_, curr_amen_area_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_publics_node_u4td3ir3jf__makeAmen_($p, curr_amen_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let curr_resi_ = [];

      let curr_resi_area_ = 0;

      if (ifn.len(parts_side_remainder_)) {

        curr_resi_area_ = 0;

        curr_resi_area_ = await exec_publics_node_u4td3ir3jf__addPartsUntilReqArea_($p, curr_resi_, parts_side_remainder_, curr_resi_area_, resi_req_area_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        neighbours_ = mfn.query.Neighbor('pg', curr_resi_);

        let resi_corners_ = await exec_publics_node_u4td3ir3jf__filterKeepCorners_($p, neighbours_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        curr_resi_area_ = await exec_publics_node_u4td3ir3jf__addParts_($p, curr_resi_, resi_corners_, curr_resi_area_);
        if ($p.terminated) {
          return mfn.getModel();
        }
      }

      let parts_remainder_ = await exec_publics_node_u4td3ir3jf__filterRem_($p, parts_, curr_resi_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      parts_remainder_ = await exec_publics_node_u4td3ir3jf__filterRem_($p, parts_remainder_, curr_amen_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_publics_node_u4td3ir3jf__makeOpen_($p, parts_remainder_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let areas_ = mfn.calc.Area(parts_remainder_);

      let curr_open_area_ = ifn.sum(areas_);

      let result_open_ = await exec_publics_node_u4td3ir3jf__printStr_($p, curr_open_area_, open_req_area_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let result_amen_ = await exec_publics_node_u4td3ir3jf__printStr_($p, curr_amen_area_, amen_req_area_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let result_resi_ = await exec_publics_node_u4td3ir3jf__printStr_($p, curr_resi_area_, resi_req_area_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      return [curr_open_area_, curr_amen_area_];
    }


    async function exec_publics_node_u4td3ir3jf__printStr_($p, curr_, req_) {

      let str_ = "current: " + ifn.numToStr(ifn.round(curr_)) + "m2,  req: " + ifn.round(req_) + " m2";

      return str_;
    }


    async function exec_publics_node_u4td3ir3jf__sort_open_($p, parts_, site_cen_, main_park_hull_) {

      let result_ = null;

      if (main_park_hull_ != null) {

        result_ = await exec_publics_node_u4td3ir3jf__sortPartsDistEdgeToPark_($p, parts_, main_park_hull_);
        if ($p.terminated) {
          return mfn.getModel();
        }
      } else {

        result_ = await exec_publics_node_u4td3ir3jf__sortPartsDistEdgeToCen_($p, parts_, site_cen_);
        if ($p.terminated) {
          return mfn.getModel();
        }
      }

      let newresult_ = [];

      for (let r_ of result_) {

        if (mfn.attrib.Get(r_[pythonList(0, r_.length)], 'type') != 'concave_corner') {

          mfn.list.Add(newresult_, r_, 'to_end');
        }
      }

      for (let r_ of result_) {

        if (mfn.attrib.Get(r_[pythonList(0, r_.length)], 'type') == 'concave_corner') {

          mfn.list.Add(newresult_, r_, 'to_end');
        }
      }

      return newresult_;
    }


    async function exec_publics_node_u4td3ir3jf__sort_amen_($p, parts_, site_cen_, main_park_hull_) {

      let result_ = null;

      if (main_park_hull_ != null) {

        result_ = await exec_publics_node_u4td3ir3jf__sortPartsDistEdgeToPark_($p, parts_, main_park_hull_);
        if ($p.terminated) {
          return mfn.getModel();
        }
      } else {

        result_ = await exec_publics_node_u4td3ir3jf__sortPartsDistEdgeToCen_($p, parts_, site_cen_);
        if ($p.terminated) {
          return mfn.getModel();
        }
      }

      return result_;
    }


    async function exec_publics_node_u4td3ir3jf__sortPartsStrightDistToCen_($p, parts_, site_cen_) {

      if (ifn.len(parts_) == 0) {

        return [];
      }

      let cens_ = mfn.calc.Centroid(parts_, 'ps_average');

      let dists_ = ifn.vecLen(ifn.vecFromTo(site_cen_, cens_));

      let sorted_ = ifn.listZip(ifn.listSort(ifn.listZip([parts_, dists_]), dists_))[pythonList(0, ifn.listZip(ifn.listSort(ifn.listZip([parts_, dists_]), dists_)).length)];

      return sorted_;
    }


    async function exec_publics_node_u4td3ir3jf__sortPartsDistEdgeToCen_($p, parts_, site_cen_) {

      if (ifn.len(parts_) == 0) {

        return [];
      }

      let unsorted_ = [];

      let sort_vals_ = [];

      for (let part_ of parts_) {

        let edge_dist_ = await exec_publics_node_u4td3ir3jf__distEdgeToCen_($p, part_, site_cen_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.list.Add(unsorted_, [part_, edge_dist_[pythonList(0, edge_dist_.length)], edge_dist_[pythonList(1, edge_dist_.length)]], 'to_end');

        mfn.list.Add(sort_vals_, edge_dist_[pythonList(1, edge_dist_.length)], 'to_end');
      }

      let sorted_ = ifn.listSort(unsorted_, sort_vals_);

      return sorted_;
    }


    async function exec_publics_node_u4td3ir3jf__sortPartsDistEdgeToPark_($p, parts_, main_park_hull_) {

      if (ifn.len(parts_) == 0) {

        return [];
      }

      let unsorted_ = [];

      let sort_vals_ = [];

      for (let part_ of parts_) {

        let edge_dist_ = await exec_publics_node_u4td3ir3jf__distEdgeToPark_($p, part_, main_park_hull_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.list.Add(unsorted_, [part_, edge_dist_[pythonList(0, edge_dist_.length)], edge_dist_[pythonList(1, edge_dist_.length)]], 'to_end');

        mfn.list.Add(sort_vals_, edge_dist_[pythonList(1, edge_dist_.length)], 'to_end');
      }

      let sorted_ = ifn.listSort(unsorted_, sort_vals_);

      sorted_ = await exec_publics_node_u4td3ir3jf__sortPartsThroughLine_($p, sorted_, main_park_hull_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      return sorted_;
    }


    async function exec_publics_node_u4td3ir3jf__sortPartsThroughLine_($p, sorted_, main_park_hull_) {

      let reorganized_ = [];

      let i_ = 0;

      let check_dict_ = {};

      while (i_ < ifn.len(sorted_)) {

        if (check_dict_[pythonList(sorted_[pythonList(i_, sorted_.length)][pythonList(0, sorted_[pythonList(i_, sorted_.length)].length)], check_dict_.length)]) {

          i_ = i_ + 1;

          continue;
        }

        check_dict_[pythonList(sorted_[pythonList(i_, sorted_.length)][pythonList(0, sorted_[pythonList(i_, sorted_.length)].length)], check_dict_.length)] = true;

        mfn.list.Add(reorganized_, sorted_[pythonList(i_, sorted_.length)], 'to_end');

        let a_ = await exec_publics_node_u4td3ir3jf__sortThroughLineComponent_($p, reorganized_, sorted_, check_dict_, i_, main_park_hull_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        i_ = i_ + 1;
      }

      return reorganized_;
    }


    async function exec_publics_node_u4td3ir3jf__sortThroughLineComponent_($p, reorganized_, sorted_, check_dict_, i_, main_park_hull_) {

      for (let j_ of ifn.range(i_ + 1, ifn.len(sorted_))) {

        if (check_dict_[pythonList(sorted_[pythonList(j_, sorted_.length)][pythonList(0, sorted_[pythonList(j_, sorted_.length)].length)], check_dict_.length)]) {

          continue;
        }

        let dists_ = mfn.calc.Distance(mfn.query.Get('ps', sorted_[pythonList(i_, sorted_.length)][pythonList(0, sorted_[pythonList(i_, sorted_.length)].length)]), mfn.query.Get('_e', sorted_[pythonList(j_, sorted_.length)][pythonList(0, sorted_[pythonList(j_, sorted_.length)].length)]), 'ps_to_e_distance');

        let checkcount_ = 0;

        for (let d_ of dists_) {

          if (d_ < 1) {

            checkcount_ = checkcount_ + 1;
          }
        }

        let mindist_ = ifn.min(dists_);

        if (checkcount_ > 1) {

          if (mfn.attrib.Get(sorted_[pythonList(j_, sorted_.length)][pythonList(0, sorted_[pythonList(j_, sorted_.length)].length)], 'type') == 'loc' && mfn.attrib.Get(sorted_[pythonList(j_, sorted_.length)][pythonList(0, sorted_[pythonList(j_, sorted_.length)].length)], 'type') == 'loc') {

            check_dict_[pythonList(sorted_[pythonList(j_, sorted_.length)][pythonList(0, sorted_[pythonList(j_, sorted_.length)].length)], check_dict_.length)] = true;

            continue;
          }

          let sortedj_ = sorted_[pythonList(j_, sorted_.length)];

          let sortedjtype_ = mfn.attrib.Get(sorted_[pythonList(j_, sorted_.length)][pythonList(0, sorted_[pythonList(j_, sorted_.length)].length)], 'type');

          let max_edge_;
          let max_dist_;
          [max_edge_, max_dist_] = await exec_publics_node_u4td3ir3jf__maxdistEdgeToPark_($p, sorted_[pythonList(j_, sorted_.length)][pythonList(0, sorted_[pythonList(j_, sorted_.length)].length)], main_park_hull_);

          let new_entry_ = [sorted_[pythonList(j_, sorted_.length)][pythonList(0, sorted_[pythonList(j_, sorted_.length)].length)], max_edge_, sorted_[pythonList(j_, sorted_.length)][pythonList(2, sorted_[pythonList(j_, sorted_.length)].length)]];

          mfn.list.Add(reorganized_, new_entry_, 'to_end');

          check_dict_[pythonList(sorted_[pythonList(j_, sorted_.length)][pythonList(0, sorted_[pythonList(j_, sorted_.length)].length)], check_dict_.length)] = true;

          if (mfn.attrib.Get(sorted_[pythonList(j_, sorted_.length)][pythonList(0, sorted_[pythonList(j_, sorted_.length)].length)], 'type') == 'off_grid0') {

            let check_ = await exec_publics_node_u4td3ir3jf__sortThroughLineComponent_($p, reorganized_, sorted_, check_dict_, j_, main_park_hull_);
            if ($p.terminated) {
              return mfn.getModel();
            }

            if (check_) {

              new_entry_[pythonList(1, new_entry_.length)] = null;
            }
          }

          return true;
        }
      }

      return false;
    }


    async function exec_publics_node_u4td3ir3jf__distEdgeToPark_($p, part_, park_) {

      let min_dist_ = Infinity;

      let min_edge_ = null;

      let acc_edges_ = mfn.query.Filter(mfn.query.Get('_e', part_), ['road', null], '==', "road_loc");

      if (ifn.len(acc_edges_) == 0) {

        acc_edges_ = mfn.query.Filter(mfn.query.Get('_e', part_), ['road', null], '==', "loc");
      }

      if (ifn.len(acc_edges_) == 0) {

        acc_edges_ = mfn.query.Filter(mfn.query.Get('_e', part_), ['road', null], '==', "sec");
      }

      if (ifn.len(acc_edges_) == 0) {

        acc_edges_ = mfn.query.Get('_e', part_);
      }

      for (let edge_ of acc_edges_) {

        let cen_ = mfn.calc.Centroid(edge_, 'ps_average');

        let posi_ = mfn.make.Position(cen_);

        let dist_ = mfn.calc.Distance(posi_, mfn.query.Get('_w', park_), 'ps_to_w_distance');

        if (dist_ < min_dist_) {

          min_dist_ = dist_;

          min_edge_ = edge_;
        }
      }

      return [min_edge_, min_dist_];
    }


    async function exec_publics_node_u4td3ir3jf__maxdistEdgeToPark_($p, part_, park_) {

      let max_dist_ = -Infinity;

      let max_edge_ = null;

      let acc_edges_ = mfn.query.Filter(mfn.query.Get('_e', part_), ['road', null], '==', "road_loc");

      if (ifn.len(acc_edges_) == 0) {

        acc_edges_ = mfn.query.Filter(mfn.query.Get('_e', part_), ['road', null], '==', "loc");
      }

      if (ifn.len(acc_edges_) == 0) {

        acc_edges_ = mfn.query.Filter(mfn.query.Get('_e', part_), ['road', null], '==', "sec");
      }

      if (ifn.len(acc_edges_) == 0) {

        acc_edges_ = mfn.query.Get('_e', part_);
      }

      for (let edge_ of acc_edges_) {

        let cen_ = mfn.calc.Centroid(edge_, 'ps_average');

        let posi_ = mfn.make.Position(cen_);

        let dist_ = mfn.calc.Distance(posi_, mfn.query.Get('_w', park_), 'ps_to_w_distance');

        if (dist_ > max_dist_) {

          max_dist_ = dist_;

          max_edge_ = edge_;
        }
      }

      return [max_edge_, max_dist_];
    }


    async function exec_publics_node_u4td3ir3jf__distEdgeToCen_($p, part_, site_cen_) {

      let min_dist_ = Infinity;

      let min_edge_ = null;

      let acc_edges_ = mfn.query.Filter(mfn.query.Get('_e', part_), ['road', null], '==', "road_loc");

      if (ifn.len(acc_edges_) == 0) {

        acc_edges_ = mfn.query.Filter(mfn.query.Get('_e', part_), ['road', null], '==', "loc");
      }

      if (ifn.len(acc_edges_) == 0) {

        acc_edges_ = mfn.query.Filter(mfn.query.Get('_e', part_), ['road', null], '==', "sec");
      }

      if (ifn.len(acc_edges_) == 0) {

        acc_edges_ = mfn.query.Get('_e', part_);
      }

      for (let edge_ of acc_edges_) {

        let cen_ = mfn.calc.Centroid(edge_, 'ps_average');

        let dist_ = ifn.distance(site_cen_, cen_);

        if (dist_ < min_dist_) {

          min_dist_ = dist_;

          min_edge_ = edge_;
        }
      }

      return [min_edge_, min_dist_];
    }


    async function exec_publics_node_u4td3ir3jf__divParts_($p, parts_) {

      let on_grid_parts_ = [];

      let on_grid_corners_ = [];

      let off_grid_parts_ = [];

      for (let part_ of parts_) {

        let type_ = mfn.attrib.Get(part_, 'type');

        if (type_ == "off_grid0") {

          mfn.list.Add(off_grid_parts_, part_, 'to_end');
        } else {
          if (type_ == "loc" || type_ == "sec") {

            mfn.list.Add(on_grid_parts_, part_, 'to_end');
          } else {

            mfn.list.Add(on_grid_corners_, part_, 'to_end');
          }
        }
      }

      return [on_grid_corners_, on_grid_parts_, off_grid_parts_];
    }


    async function exec_publics_node_u4td3ir3jf__reorderParts_($p, parts_) {

      let result_ = await exec_publics_node_u4td3ir3jf__divParts_($p, parts_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let reordered_parts_ = [];

      mfn.list.Add(reordered_parts_, result_[pythonList(0, result_.length)][pythonList(0, result_[pythonList(0, result_.length)].length)], 'to_end');

      reordered_parts_ = ifn.listJoin(reordered_parts_, result_[pythonList(1, result_.length)].slice(0, 2));

      reordered_parts_ = ifn.listJoin(reordered_parts_, result_[pythonList(2, result_.length)]);

      reordered_parts_ = ifn.listJoin(reordered_parts_, result_[pythonList(0, result_.length)].slice(1, -1));

      reordered_parts_ = ifn.listJoin(reordered_parts_, result_[pythonList(1, result_.length)].slice(2));

      mfn.list.Add(reordered_parts_, result_[pythonList(0, result_.length)][pythonList(-1, result_[pythonList(0, result_.length)].length)], 'to_end');

      return reordered_parts_;

      return parts_;

      let a_ = 1;
    }


    async function exec_publics_node_u4td3ir3jf__reorderParts1_($p, parts_) {

      if (ifn.len(parts_) < 1) {

        return [];
      }

      let og_part_list_ = [];

      let other_part_list_ = [];

      for (let part_ of parts_) {

        let part_type_ = mfn.attrib.Get(part_, 'type');

        let partcoord_ = mfn.attrib.Get(mfn.query.Get('ps', part_), 'xyz');

        let part_midpt_ = ifn.vecDiv(ifn.vecSum(partcoord_), ifn.len(partcoord_));

        if (part_type_ == 'off_grid' || part_type_ == 'off_grid0') {

          mfn.list.Add(og_part_list_, [part_, 'off_grid', part_midpt_, []], 'to_end');
        } else {

          mfn.list.Add(other_part_list_, [part_, part_type_, part_midpt_], 'to_end');
        }
      }

      for (let part_item_ of other_part_list_) {

        let shortest_item_ = [1000000000, null];

        for (let og_part_item_ of og_part_list_) {

          let dist_vec_ = ifn.vecSub(part_item_[pythonList(2, part_item_.length)], og_part_item_[pythonList(2, og_part_item_.length)]);

          let dist2_ = ifn.vecDot(dist_vec_, dist_vec_);

          if (shortest_item_[pythonList(0, shortest_item_.length)] > dist2_) {

            shortest_item_ = [dist2_, og_part_item_];
          }
        }

        if (part_item_[pythonList(1, part_item_.length)] == 'loc') {

          mfn.list.Add(shortest_item_[pythonList(1, shortest_item_.length)][pythonList(3, shortest_item_[pythonList(1, shortest_item_.length)].length)], part_item_[pythonList(0, part_item_.length)], 'to_start');
        } else {

          mfn.list.Add(shortest_item_[pythonList(1, shortest_item_.length)][pythonList(3, shortest_item_[pythonList(1, shortest_item_.length)].length)], part_item_[pythonList(0, part_item_.length)], 'to_end');
        }
      }

      return og_part_list_;
    }


    async function exec_publics_node_u4td3ir3jf__reorderParts1_copy_($p, parts_) {

      if (ifn.len(parts_) < 1) {

        return [[], []];
      }

      let x_ = "==============================";
      printFunc($p.console, `x`, x_);

      let edgex_ = mfn.query.Get('_e', parts_[pythonList(0, parts_.length)])[pythonList(0, mfn.query.Get('_e', parts_[pythonList(0, parts_.length)]).length)];
      printFunc($p.console, `edgex`, edgex_);

      let edgey_ = mfn.query.Get('_e', parts_[pythonList(0, parts_.length)])[pythonList(1, mfn.query.Get('_e', parts_[pythonList(0, parts_.length)]).length)];
      printFunc($p.console, `edgey`, edgey_);

      let rayx_ = mfn.calc.Ray(edgex_);
      printFunc($p.console, `rayx`, rayx_);

      let rayy_ = mfn.calc.Ray(edgey_);
      printFunc($p.console, `rayy`, rayy_);

      let coefflist_ = [];

      let partlist_ = [];

      let cx_ = 0;

      if (ifn.abs(rayx_[pythonList(1, rayx_.length)][pythonList(0, rayx_[pythonList(1, rayx_.length)].length)]) < ifn.abs(rayx_[pythonList(1, rayx_.length)][pythonList(1, rayx_[pythonList(1, rayx_.length)].length)])) {

        cx_ = 1;
      }

      let cy_ = 0;

      if (ifn.abs(rayy_[pythonList(1, rayy_.length)][pythonList(0, rayy_[pythonList(1, rayy_.length)].length)]) < ifn.abs(rayy_[pythonList(1, rayy_.length)][pythonList(1, rayy_[pythonList(1, rayy_.length)].length)])) {

        cy_ = 1;
      }

      for (let part_ of parts_) {

        let partcoord_ = mfn.attrib.Get(mfn.query.Get('ps', part_), 'xyz');
        printFunc($p.console, `partcoord`, partcoord_);

        let midpt_ = ifn.vecDiv(ifn.vecSum(partcoord_), ifn.len(partcoord_));
        printFunc($p.console, `midpt`, midpt_);

        let part_type_ = mfn.attrib.Get(part_, 'type');
      }

      return [];
    }


    async function exec_publics_node_u4td3ir3jf__addPartsUntilReqArea_($p, curr_list_, parts_to_add_, curr_area_, req_area_) {

      req_area_ = req_area_;

      for (let part_ of parts_to_add_) {

        if (curr_area_ < req_area_) {

          mfn.list.Add(curr_list_, part_, 'to_end');

          let area_ = mfn.calc.Area(part_);

          curr_area_ = curr_area_ + area_;
        } else {

          break;
        }
      }

      return curr_area_;
    }


    async function exec_publics_node_u4td3ir3jf__addParts_($p, curr_list_, parts_to_add_, curr_area_) {

      for (let part_ of parts_to_add_) {

        let area_ = mfn.calc.Area(part_);

        mfn.list.Add(curr_list_, part_, 'to_end');

        curr_area_ = curr_area_ + area_;
      }

      return curr_area_;
    }


    async function exec_publics_node_u4td3ir3jf__makeOpen_($p, parts_) {

      for (let part_ of parts_) {

        let type_ = mfn.attrib.Get(part_, 'type');

        if (type_ == 'concave_corner') {

          type_ = 'off_grid2';
        }

        if (!ifn.strEnds(type_, "_os")) {

          mfn.attrib.Set(part_, `type`, type_ + "_os");

          mfn.attrib.Set(part_, `class`, "open_main");
        }

        mfn.visualize.Color(part_, [0, 0.5, 0]);
      }
    }


    async function exec_publics_node_u4td3ir3jf__makeAmen_($p, parts_) {

      for (let part_ of parts_) {

        let type_ = mfn.attrib.Get(part_, 'type');

        if (type_ == 'concave_corner') {

          type_ = 'off_grid2';
        }

        if (!ifn.strEnds(type_, "_am")) {

          mfn.attrib.Set(part_, `type`, type_ + "_am");
        }

        mfn.visualize.Color(part_, [0.5, 0, 0]);
      }
    }


    async function exec_publics_node_u4td3ir3jf__filterKeepCorners_($p, parts_) {

      let corners_ = [];

      for (let part_ of parts_) {

        if (mfn.attrib.Get(part_, 'type') == "loc_loc" || mfn.attrib.Get(part_, 'type') == "sec_loc") {

          mfn.list.Add(corners_, part_, 'to_end');
        }
      }

      return corners_;
    }


    async function exec_publics_node_u4td3ir3jf__filterKeepOg_($p, parts_) {

      let og_ = [];

      for (let part_ of parts_) {

        let type_ = mfn.attrib.Get(part_, 'type');

        if (type_ == "off_grid0" || type_ == "off_grid1") {

          mfn.list.Add(og_, part_, 'to_end');
        }

        if (type_ == "concave_corner") {

          mfn.list.Add(og_, part_, 'to_start');
        }
      }

      return og_;
    }


    async function exec_publics_node_u4td3ir3jf__filterRem_($p, parts_, parts_to_rem_) {

      let filtered_ = [];

      for (let part_ of parts_) {

        if (!ifn.listHas(parts_to_rem_, part_)) {

          mfn.list.Add(filtered_, part_, 'to_end');
        }
      }

      return filtered_;
    }


    async function exec_publics_node_u4td3ir3jf__filterThroughpart_($p, partgroups_) {

      let og_ = [];

      let startcheck_ = true;

      let endcheck_ = true;

      let endpart_ = null;

      for (let parts_ of partgroups_) {

        endcheck_ = true;

        for (let part_ of parts_) {

          if (mfn.attrib.Get(part_, 'type') == 'loc') {

            if (startcheck_) {

              mfn.list.Add(og_, part_, 'to_start');

              startcheck_ = false;
            } else {
              if (endcheck_) {

                endpart_ = part_;
              }
            }
          } else {
            if (mfn.attrib.Get(part_, 'type') == 'off_grid0') {

              mfn.list.Add(og_, part_, 'to_end');
            } else {
              if (mfn.attrib.Get(part_, 'type') == 'loc_loc') {

                endcheck_ = false;
              }
            }
          }
        }
      }

      return og_;
    }


    async function exec_publics_node_u4td3ir3jf__coptAttribs_($p, from_part_, to_part_) {

      mfn.attrib.Set(to_part_, `type`, mfn.attrib.Get(from_part_, 'type'));

      mfn.attrib.Set(to_part_, `site`, mfn.attrib.Get(from_part_, 'site'));

      mfn.attrib.Set(to_part_, `block_type`, mfn.attrib.Get(from_part_, 'block_type'));

      mfn.attrib.Set(to_part_, `block_id`, mfn.attrib.Get(from_part_, 'block_id'));

      mfn.attrib.Set(to_part_, `class`, mfn.attrib.Get(from_part_, 'class'));
    }


    async function exec_publics_node_u4td3ir3jf__transferEdgeAttribs_($p, from_edges_, to_edges_) {

      for (let to_edge_ of to_edges_) {

        let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

        let from_edge_ = await exec_publics_node_u4td3ir3jf__touchingEdge_($p, from_edges_, cen_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (from_edge_ != null) {

          let val_ = mfn.attrib.Get(from_edge_, 'road');

          if (val_ != undefined) {

            mfn.attrib.Set(to_edge_, `road`, val_);
          }
        }
      }
    }


    async function exec_publics_node_u4td3ir3jf__transferEdgeAttribsBtwTouchingParts_($p, parts_) {

      let edges_ = mfn.query.Filter(mfn.query.Get('_e', parts_), ['road', null], '==', null);

      for (let to_edge_ of edges_) {

        if (mfn.attrib.Get(to_edge_, 'road') == undefined) {

          let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

          let idx_ = ifn.listFind(edges_, to_edge_);

          let from_edges_ = ifn.listJoin(edges_.slice(0, idx_), edges_.slice(idx_ + 1));

          let from_edge_ = await exec_publics_node_u4td3ir3jf__touchingEdge_($p, from_edges_, cen_);
          if ($p.terminated) {
            return mfn.getModel();
          }

          if (from_edge_ != null) {

            let part_type_ = mfn.attrib.Get(mfn.query.Get('pg', from_edge_)[pythonList(0, mfn.query.Get('pg', from_edge_).length)], 'type');

            if (part_type_ != undefined) {

              mfn.attrib.Set(to_edge_, `road`, part_type_);
            }
          }
        }
      }
    }


    async function exec_publics_node_u4td3ir3jf__touchingEdge_($p, edges_, xyz_) {

      for (let edge_ of edges_) {

        let r_ = mfn.calc.Ray(edge_);

        let d_ = ifn.distance(xyz_, ifn.project(xyz_, r_, 0));

        if (d_ < 0.1) {

          return edge_;
        }
      }

      return null;
    }


    async function exec_publics_node_u4td3ir3jf__distEdgeToPark_copy_($p, part_, park_) {

      let min_dist_ = Infinity;

      let min_edge_ = null;

      let acc_edges_ = mfn.query.Filter(mfn.query.Get('_e', part_), ['road', null], '==', "road_loc");

      if (ifn.len(acc_edges_) == 0) {

        acc_edges_ = mfn.query.Filter(mfn.query.Get('_e', part_), ['road', null], '==', "loc");
      }

      if (ifn.len(acc_edges_) == 0) {

        acc_edges_ = mfn.query.Filter(mfn.query.Get('_e', part_), ['road', null], '==', "sec");
      }

      if (ifn.len(acc_edges_) == 0) {

        acc_edges_ = mfn.query.Get('_e', part_);
      }

      for (let edge_ of acc_edges_) {

        let cen_ = mfn.calc.Centroid(edge_, 'ps_average');

        let posi_ = mfn.make.Position(cen_);

        let dist_ = mfn.calc.Distance(posi_, mfn.query.Get('_w', park_), 'ps_to_w_distance');

        if (dist_ < min_dist_) {

          min_dist_ = dist_;

          min_edge_ = edge_;
        }
      }

      return [min_edge_, min_dist_];
    }

    async function exec_publics_node_u4td3ir3jf($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: open and amen', '__null__')
      }


      let att_check_ = false;

      let pg_attr_ = mfn.attrib.Discover('pg');

      for (let attr_ of pg_attr_) {

        if (attr_['name'] == 'cluster_id') {

          att_check_ = true;
        }
      }

      if (!att_check_) {

        mfn.attrib.Add('pg', 'string', 'cluster_id');
      }

      for (let site_name_ of mfn.attrib.Get(null, 'site_names')) {

        let oldparts_ = mfn.query.Filter(mfn.query.Get('pg', null), ['site', null], '==', site_name_);

        let site_cen_ = mfn.attrib.Get(null, site_name_ + "_cen");

        let block_id_ = mfn.attrib.Get(null, site_name_ + "_cen_block_id");

        let parts_ = mfn.make.Clone(oldparts_);

        parts_ = mfn.query.Get('pg', parts_);

        let site_area_ = mfn.attrib.Get(null, site_name_ + "_area");

        let req_open_area_ = site_area_ * (OPEN_PERCENT_ / 100);
        printFunc($p.console, `req_open_area`, req_open_area_);

        let req_amen_area_ = site_area_ * (AMEN_PERCENT_ / 100);

        let result_ = await exec_publics_node_u4td3ir3jf_getCenBlk_($p, parts_, block_id_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        let parts_cen_ = result_[pythonList(0, result_.length)];

        let parts_other_ = result_[pythonList(1, result_.length)];

        result_ = await exec_publics_node_u4td3ir3jf_processCenBlk_($p, parts_cen_, site_cen_, req_open_area_, req_amen_area_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        let curr_open_area_ = 0;

        let curr_amen_area_ = 0;

        let main_park_hull_ = null;

        if (result_ != null) {

          curr_open_area_ = result_[pythonList(0, result_.length)];

          curr_amen_area_ = result_[pythonList(1, result_.length)];

          main_park_hull_ = result_[pythonList(2, result_.length)];
        }

        let loc_parts_ = mfn.query.Filter(parts_other_, ['type', null], '==', "loc");

        let loc_loc_parts_ = mfn.query.Filter(parts_other_, ['type', null], '==', "loc_loc");

        let sec_parts_ = mfn.query.Filter(parts_other_, ['type', null], '==', "sec");

        parts_other_ = mfn.query.Get('pg', parts_other_);

        let ls_parts_ = ifn.listJoin(loc_parts_, loc_loc_parts_, sec_parts_);

        let og_parts_ = await exec_publics_node_u4td3ir3jf__filterKeepOg_($p, parts_other_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        let og_parts_edges_ = await exec_publics_node_u4td3ir3jf__sort_open_($p, og_parts_, site_cen_, main_park_hull_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        og_parts_edges_ = ifn.listRev(og_parts_edges_);

        result_ = await exec_publics_node_u4td3ir3jf_insertIntoBlks_($p, og_parts_edges_, ls_parts_, curr_open_area_, req_open_area_, true);
        if ($p.terminated) {
          return mfn.getModel();
        }

        let added_parts_ = result_[pythonList(0, result_.length)];

        let new_loc_parts_ = result_[pythonList(1, result_.length)];

        curr_open_area_ = result_[pythonList(2, result_.length)];

        await exec_publics_node_u4td3ir3jf__makeOpen_($p, added_parts_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        parts_other_ = mfn.query.Get('pg', parts_other_);

        ls_parts_ = ifn.listJoin(mfn.query.Get('pg', ls_parts_), new_loc_parts_);

        og_parts_ = await exec_publics_node_u4td3ir3jf__filterKeepOg_($p, parts_other_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        og_parts_ = await exec_publics_node_u4td3ir3jf__filterRem_($p, og_parts_, added_parts_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        og_parts_edges_ = await exec_publics_node_u4td3ir3jf__sort_amen_($p, og_parts_, site_cen_, main_park_hull_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        result_ = await exec_publics_node_u4td3ir3jf_insertIntoBlks_($p, og_parts_edges_, ls_parts_, curr_amen_area_, req_amen_area_, true);
        if ($p.terminated) {
          return mfn.getModel();
        }

        added_parts_ = result_[pythonList(0, result_.length)];

        let new_ls_parts_ = result_[pythonList(1, result_.length)];

        curr_amen_area_ = result_[pythonList(2, result_.length)];

        await exec_publics_node_u4td3ir3jf__makeAmen_($p, added_parts_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (main_park_hull_ != null) {

          mfn.edit.Delete(main_park_hull_, 'delete_selected');
        }
      }
    }


    async function exec_publics_node_qfb5rsev7vl_applyColours_($p, colors_dict_, pgons_) {

      let pgons_lists_ = {};

      let cols_ = {};

      let keys_ = [];

      for (let pgon_ of pgons_) {

        let type_ = mfn.attrib.Get(pgon_, 'type');

        let col_ = colors_dict_[pythonList(type_, colors_dict_.length)];

        if (col_ != undefined) {

          let key_ = ifn.string(col_);

          let pgons_list_ = pgons_lists_[pythonList(key_, pgons_lists_.length)];

          if (pgons_list_ == undefined) {

            pgons_list_ = [];

            pgons_lists_[pythonList(key_, pgons_lists_.length)] = pgons_list_;

            cols_[pythonList(key_, cols_.length)] = col_;

            mfn.list.Add(keys_, key_, 'to_end');
          }

          mfn.list.Add(pgons_list_, pgon_, 'to_end');
        }
      }

      for (let key_ of keys_) {

        pgons_ = pgons_lists_[pythonList(key_, pgons_lists_.length)];

        let col_ = cols_[pythonList(key_, cols_.length)];

        mfn.visualize.Color(pgons_, col_);
      }
    }

    async function exec_publics_node_qfb5rsev7vl($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: colors', '__null__')
      }


      let colors_dict_bldg_ = {
        "art_art": [0.8, 0.51, 0.48],
        "art_sec": [0.85, 0.62, 0.6],
        "art_loc": [0.9, 0.73, 0.72],
        "art": [0.92, 0.79, 0.78],
        "sec_sec": [0.98, 0.81, 0.69],
        "sec_loc": [0.98, 0.86, 0.78],
        "sec": [0.98, 0.91, 0.86],
        "loc_loc": [0.98, 0.92, 0.73],
        "loc": [0.98, 0.95, 0.83],
        "off_grid": [0.81, 0.9, 0.95],
        "off_grid0": [0.81, 0.9, 0.95],
        "off_grid1": [0.81, 0.9, 0.95],
        "off_grid2": [0.68, 0.83, 0.9],
        "am": [0.85, 0.6, 0.66],
        "sec_am": [0.85, 0.6, 0.66],
        "sec_sec_am": [0.85, 0.6, 0.66],
        "sec_loc_am": [0.85, 0.6, 0.66],
        "loc_am": [0.85, 0.6, 0.66],
        "loc_loc_am": [0.85, 0.6, 0.66],
        "off_grid0_am": [0.85, 0.6, 0.66],
        "off_grid1_am": [0.85, 0.6, 0.66],
        "off_grid2_am": [0.85, 0.6, 0.66],
        "green0": [0.86, 0.95, 0.81],
        "green1": [0.86, 0.95, 0.81],
        "green2": [0.86, 0.95, 0.81],
        "os": [0.75, 0.9, 0.67],
        "sec_os": [0.75, 0.9, 0.67],
        "sec_sec_os": [0.75, 0.9, 0.67],
        "sec_loc_os": [0.75, 0.9, 0.67],
        "loc_os": [0.75, 0.9, 0.67],
        "loc_loc_os": [0.75, 0.9, 0.67],
        "off_grid0_os": [0.75, 0.9, 0.67],
        "off_grid1_os": [0.75, 0.9, 0.67],
        "off_grid2_os": [0.75, 0.9, 0.67],
        "corner_park": [0.75, 0.9, 0.67],
        "internal path, any cluster": [0.93, 1, 0.9],
        "path0": [0.93, 1, 0.9],
        "path1": [0.93, 1, 0.9],
        "path2": [0.93, 1, 0.9],
        "Cluster, access path": [1, 0.98, 0.92],
        "entr0": [1, 0.98, 0.92],
        "entr1": [1, 0.98, 0.92],
        "og cluster, 2nd": [0.68, 0.83, 0.9],
        "Public streets": [1, 1, 1]
      };
      printFunc($p.console, `colors_dict_bldg`, colors_dict_bldg_);

      colors_dict_bldg_["off_grid_too_small"] = [0.81, 0.9, 0.95];

      colors_dict_bldg_["too_deep"] = [0.81, 0.9, 0.95];

      colors_dict_bldg_["concave_corner"] = [0.81, 0.9, 0.95];

      colors_dict_bldg_["leftover"] = [0.81, 0.9, 0.95];

      await exec_publics_node_qfb5rsev7vl_applyColours_($p, colors_dict_bldg_, mfn.query.Get('pg', null));
      if ($p.terminated) {
        return mfn.getModel();
      }
    }


    async function exec_publics_node_b80cwerw0s($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: set input', '__null__')
      }


      mfn.attrib.Set(null, `neighbourhood_input`, {
        "road_art_w": ROAD_ART_W_,
        "road_sec_w": ROAD_SEC_W_,
        "road_loc_w": ROAD_LOC_W_,
        "part_art_d": PART_ART_D_,
        "part_sec_d": PART_SEC_D_,
        "part_loc_d": PART_LOC_D_,
        "part_og_d": PART_OG_D_,
        "part_og_w": PART_OG_W_,
        "plot_art_w": PLOT_ART_W_,
        "plot_sec_w": PLOT_SEC_W_,
        "plot_loc_w": PLOT_LOC_W_,
        "blk_art_num_og_d": BLK_ART_NUM_OG_D_,
        "blk_art_num_og_w": BLK_ART_NUM_OG_W_,
        "blk_sec_num_og_d": BLK_SEC_NUM_OG_D_,
        "blk_sec_num_og_w": BLK_SEC_NUM_OG_W_,
        "blk_loc_num_og_d": BLK_LOC_NUM_OG_D_,
        "blk_loc_num_og_w": BLK_LOC_NUM_OG_W_,
        "path_w": PATH_W_,
        "open_percent": OPEN_PERCENT_,
        "amen_percent": AMEN_PERCENT_,
        "pavement_w": PAVEMENT_W_,
        "add_trees": ADD_TREES_,
        "tree_spacing": TREE_SPACING_,
        "tree_height_start": TREE_HEIGHT_START_,
        "tree_height_max": TREE_HEIGHT_MAX_
      });
    }


    async function exec_publics_node_stckcsqz8ph_debugCheckOpenAmen_($p,) {

      let all_pg_ = mfn.query.Get('pg', null);

      let all_pg_types_ = mfn.attrib.Get(all_pg_, 'type');

      let open_pg_ = [];

      let amen_pg_ = [];

      for (let i_ of ifn.range(ifn.len(all_pg_))) {

        if (ifn.strEnds(all_pg_types_[pythonList(i_, all_pg_types_.length)], '_os')) {

          mfn.list.Add(open_pg_, all_pg_[pythonList(i_, all_pg_.length)], 'to_end');
        }

        if (ifn.strEnds(all_pg_types_[pythonList(i_, all_pg_types_.length)], '_am')) {

          mfn.list.Add(amen_pg_, all_pg_[pythonList(i_, all_pg_.length)], 'to_end');
        }
      }

      let open_areas_ = mfn.calc.Area(open_pg_);

      let amen_areas_ = mfn.calc.Area(amen_pg_);

      let open_area_ = ifn.sum(open_areas_);
      printFunc($p.console, `open_area`, open_area_);

      let open_percentage_ = open_area_ * 100 / mfn.attrib.Get(null, 'site_area');
      printFunc($p.console, `open_percentage`, open_percentage_);

      let amen_area_ = ifn.sum(amen_areas_);
      printFunc($p.console, `amen_area`, amen_area_);

      let amen_percentage_ = amen_area_ * 100 / mfn.attrib.Get(null, 'site_area');
      printFunc($p.console, `amen_percentage`, amen_percentage_);
    }


    async function exec_publics_node_stckcsqz8ph_debugCheckAreas_($p, data_) {

      let check_col_l_ = await exec_publics_node_stckcsqz8ph_checkTotalsColL_($p, data_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let check_col_m_ = await exec_publics_node_stckcsqz8ph_checkTotalsColM_($p, data_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let check_col_o_ = await exec_publics_node_stckcsqz8ph_checkTotalsColO_($p, data_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let check_col_g_ = await exec_publics_node_stckcsqz8ph_checkTotalsColG_($p, data_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let deubg_ = "===";
      printFunc($p.console, `deubg`, deubg_);

      let check_cols_lmo_total_ = check_col_l_ + check_col_m_ + check_col_o_;
      printFunc($p.console, `check_cols_lmo_total`, check_cols_lmo_total_);

      let check_blocks_area_ = mfn.attrib.Get(null, 'blocks_area');
      printFunc($p.console, `check_blocks_area`, check_blocks_area_);

      let check_blocks_missing_ = check_blocks_area_ - check_cols_lmo_total_;
      printFunc($p.console, `check_blocks_missing`, check_blocks_missing_);

      let check_blocks_error_ = ifn.string(100 * (check_blocks_missing_ / check_blocks_area_)) + "%";
      printFunc($p.console, `check_blocks_error`, check_blocks_error_);

      deubg_ = "===";
      printFunc($p.console, `deubg`, deubg_);

      let check_cols_lmog_total_ = check_col_l_ + check_col_m_ + check_col_o_ + check_col_g_;
      printFunc($p.console, `check_cols_lmog_total`, check_cols_lmog_total_);

      let check_site_area_ = mfn.attrib.Get(null, 'site_area');
      printFunc($p.console, `check_site_area`, check_site_area_);

      let check_site_missing_ = check_site_area_ - check_cols_lmog_total_;
      printFunc($p.console, `check_site_missing`, check_site_missing_);

      let check_site_error_ = ifn.string(100 * (check_site_missing_ / check_site_area_)) + "%";
      printFunc($p.console, `check_site_error`, check_site_error_);
    }


    async function exec_publics_node_stckcsqz8ph_checkTotalsColL_($p, data_) {

      let total_ = 0;

      total_ = total_ + data_["open_loc_loc_area"];

      total_ = total_ + data_["open_loc_area"];

      total_ = total_ + data_["open_og_clus0_on_art_area"];

      total_ = total_ + data_["open_og_clus0_on_sec_area"];

      total_ = total_ + data_["open_og_clus0_on_loc_area"];

      total_ = total_ + data_["open_og_clus1_on_loc_area"];

      total_ = total_ + data_["open_og_clus2_on_loc_area"];

      return total_;
    }


    async function exec_publics_node_stckcsqz8ph_checkTotalsColM_($p, data_) {

      let total_ = 0;

      total_ = total_ + data_["amen_loc_loc_area"];

      total_ = total_ + data_["amen_loc_area"];

      total_ = total_ + data_["amen_og_clus0_on_art_area"];

      total_ = total_ + data_["amen_og_clus0_on_sec_area"];

      total_ = total_ + data_["amen_og_clus0_on_loc_area"];

      total_ = total_ + data_["amen_og_clus1_on_loc_area"];

      total_ = total_ + data_["amen_og_clus2_on_loc_area"];

      return total_;
    }


    async function exec_publics_node_stckcsqz8ph_checkTotalsColO_($p, data_) {

      let total_ = 0;

      total_ = total_ + data_["lot_art_art_area"];

      total_ = total_ + data_["lot_art_sec_area"];

      total_ = total_ + data_["lot_art_loc_area"];

      total_ = total_ + data_["lot_art_area"];

      total_ = total_ + data_["lot_sec_sec_area"];

      total_ = total_ + data_["lot_sec_loc_area"];

      total_ = total_ + data_["lot_sec_area"];

      total_ = total_ + data_["lot_loc_loc_area"];

      total_ = total_ + data_["lot_loc_area"];

      total_ = total_ + data_["og_clus0_on_art_area"];

      total_ = total_ + data_["og_clus0_on_sec_area"];

      total_ = total_ + data_["og_clus0_on_loc_area"];

      total_ = total_ + data_["og_clus1_on_loc_area"];

      total_ = total_ + data_["og_clus2_on_loc_area"];

      return total_;
    }


    async function exec_publics_node_stckcsqz8ph_checkTotalsColG_($p, data_) {

      let total_ = 0;

      total_ = total_ + data_["road_area_art"];

      total_ = total_ + data_["road_area_sec"];

      total_ = total_ + data_["road_area_loc"];

      return total_;
    }

    async function exec_publics_node_stckcsqz8ph($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: End', '__null__')
      }


      let data_ = {};

      await exec_publics_node_stckcsqz8ph_debugCheckOpenAmen_($p);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await mfn.io.Export(null, "public.sim", 'sim', 'Save to Local Storage');

      let pgons_ = mfn.query.Get('pg', null);

      for (let pg_ of pgons_) {

        let rgb_ = mfn.attrib.Get(mfn.query.Get('_v', pg_)[pythonList(0, mfn.query.Get('_v', pg_).length)], 'rgb');

        if (rgb_) {

          let s_ = [];

          for (let val_ of rgb_) {

            mfn.list.Add(s_, ifn.round(val_ * 255), 'to_end');
          }

          mfn.attrib.Set(pg_, `color`, 'rgb(' + s_[pythonList(0, s_.length)] + ',' + s_[pythonList(1, s_.length)] + ',' + s_[pythonList(2, s_.length)] + ')');
        } else {

          mfn.attrib.Set(pg_, `color`, '#ffffff');
        }
      }

      let __return_value__ = data_;
      return __return_value__;
    }

    var merged;
    let ssid_exec_publics_node_d36n03wr8tj = mfn.model.snapshotGetActive();

    let ssid_exec_publics_node_huiuxvcjtij = ssid_exec_publics_node_d36n03wr8tj;

    await exec_publics_node_huiuxvcjtij($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_publics_node_u4td3ir3jf = ssid_exec_publics_node_huiuxvcjtij;

    await exec_publics_node_u4td3ir3jf($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_publics_node_qfb5rsev7vl = ssid_exec_publics_node_u4td3ir3jf;

    await exec_publics_node_qfb5rsev7vl($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_publics_node_b80cwerw0s = ssid_exec_publics_node_qfb5rsev7vl;

    await exec_publics_node_b80cwerw0s($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_publics_node_stckcsqz8ph = mfn.model.snapshotNext([ssid_exec_publics_node_b80cwerw0s]);

    return await exec_publics_node_stckcsqz8ph($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);
  }


  function pythonList(x, l) {
    if (x < 0) {
      return x + l;
    }
    return x;
  }

  function printFunc(_console, name, value) {
    let val;
    let padding_style = 'padding: 2px 0px 2px 10px;';
    if (!value) {
      val = value;
    } else if (value === '__null__') {
      _console.push('<p style="' + padding_style + '"><b><i>_ ' + name + '</i></b></p>');
      return value;
    } else if (typeof value === 'number' || value === undefined) {
      val = value;
    } else if (typeof value === 'string') {
      val = '"' + value.replace(/\n/g, '<br>') + '"';
    } else if (value.constructor === [].constructor) {
      let __list_check__ = false;
      let __value_strings__ = [];
      for (const __item__ of value) {
        if (!__item__) {
          __value_strings__.push('' + __item__);
          continue;
        }
        if (__item__.constructor === [].constructor || __item__.constructor === {}.constructor) {
          __list_check__ = true;
        }
        __value_strings__.push(JSON.stringify(__item__).replace(/,/g, ', '));
      }
      if (__list_check__) {
        padding_style = 'padding: 2px 0px 0px 10px;';
        val = '[<p style="padding: 0px 0px 2px 40px;">' +
          __value_strings__.join(',</p><p style="padding: 0px 0px 2px 40px;">') +
          '</p><p style="padding: 0px 0px 2px 30px;">]</p>';
      } else {
        val = '[' + __value_strings__.join(', ') + ']';
      }
    } else if (value.constructor === {}.constructor) {
      let __list_check__ = false;
      let __value_strings__ = [];
      for (const __item__ in value) {
        const __value__ = value[__item__];
        if (!__value__) {
          __value_strings__.push('\<b>"' + __item__ + '\"</b>' + ': ' + __value__);
          continue;
        }
        if (__value__.constructor === [].constructor || __value__.constructor === {}.constructor) {
          __list_check__ = true;
        }
        __value_strings__.push('\<b>"' + __item__ + '\"</b>' + ': ' + JSON.stringify(__value__).replace(/,/g, ', '));
      }
      if (__list_check__) {
        padding_style = 'padding: 2px 0px 0px 10px;';
        val = '{<p style="padding: 0px 0px 2px 40px;">' +
          __value_strings__.join(',</p><p style="padding: 0px 0px 2px 40px;">') +
          '</p><p style="padding: 0px 0px 2px 30px;">}</p>';
      } else {
        val = '{' + __value_strings__.join(', ') + '}';
      }
    } else {
      val = value;
    }
    _console.push('<p style="' + padding_style + '"><b><i>_ ' + name + '</i></b>  = ' + val + '</p>');
    return val;
  }


  const $p = {};
  if (__model__) {
    mfn.io.ImportData(__model__, 'gi');
  }
  mfn.getModel().debug = true;
  $p["console"] = [];
  $p["modules"] = mfn;
  $p["curr_ss"] = {};
  const result = await exec_publics($p, IN_MODEL, FORCE_STOP_MOBIUS, ROAD_ART_W, ROAD_SEC_W, ROAD_LOC_W, PART_ART_D, PART_SEC_D, PART_LOC_D, PART_OG_D, PART_OG_W, PLOT_ART_W, PLOT_SEC_W, PLOT_LOC_W, BLK_ART_NUM_OG_D, BLK_ART_NUM_OG_W, BLK_SEC_NUM_OG_D, BLK_SEC_NUM_OG_W, BLK_LOC_NUM_OG_D, BLK_LOC_NUM_OG_W, PATH_W, OPEN_PERCENT, AMEN_PERCENT, PAVEMENT_W, ADD_TREES, TREE_SPACING, TREE_HEIGHT_START, TREE_HEIGHT_MAX);
  if (result === mfn.getModel()) {
    return { "model": mfn.getModel(), "result": null };
  }
  return { "model": mfn.getModel(), "result": result };
  /** * **/

}

module.exports = publics;
