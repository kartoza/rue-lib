/**
 * to use this code: import cluster1 from this js file as well as the GI module
 * run cluster1 with the GI module as input along with other start node input
 * e.g.:
 * const cluster1 = require('./cluster1.js').cluster1
 * const module = require('gi-module')
 * const result = await cluster1(module, start_input_1, start_input_2, ...);
 *
 * returns: a json object:
 *   _ result.model -> gi model of the flowchart
 *   _ result.result -> returned output of the flowchart, if the flowchart does not return any value, result.result is the model of the flowchart
 */

// Parameter: {"name":"IN_MODEL","value":"street.sim","type":0}
// Parameter: {"name":"FORCE_STOP_MOBIUS","value":"[false]","type":0}
// Parameter: {"name":"ROAD_ART_W","value":"20","type":1,"min":"10","max":"30","step":"1"}
// Parameter: {"name":"ROAD_SEC_W","value":"15","type":1,"min":"10","max":"30","step":"1"}
// Parameter: {"name":"ROAD_LOC_W","value":10,"type":1,"min":"5","max":"30","step":"1"}
// Parameter: {"name":"PART_ART_D","value":"40","type":1,"min":"5","max":"50","step":"1"}
// Parameter: {"name":"PART_SEC_D","value":30,"type":1,"min":"5","max":"50","step":"1"}
// Parameter: {"name":"PART_LOC_D","value":20,"type":1,"min":"5","max":"50","step":"1"}
// Parameter: {"name":"PART_OG_D","value":"40","type":1,"min":"5","max":"50","step":"1"}
// Parameter: {"name":"PART_OG_W","value":"30","type":1,"min":"5","max":"50","step":"1"}
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
// Parameter: {"name":"OPEN_PERCENT","value":"2","type":1,"min":"0","max":"10","step":"0.1"}
// Parameter: {"name":"AMEN_PERCENT","value":"8","type":1,"min":"0","max":"15","step":"0.1"}
// Parameter: {"name":"PAVEMENT_W","value":"3","type":1,"min":"1","max":"10","step":"0.5"}
// Parameter: {"name":"ADD_TREES","value":true,"type":2}
// Parameter: {"name":"TREE_SPACING","value":12,"type":1,"min":"6","max":"30","step":"1"}
// Parameter: {"name":"TREE_HEIGHT_START","value":"8","type":1,"min":"5","max":"15","step":"1"}
// Parameter: {"name":"TREE_HEIGHT_MAX","value":"20","type":1,"min":"10","max":"30","step":"1"}


const mfn = require('@design-automation/mobius-sim-funcs').Funcs();
const ifn = require('@design-automation/mobius-inline-funcs').InlineClass(true);

async function cluster1(IN_MODEL, FORCE_STOP_MOBIUS, ROAD_ART_W, ROAD_SEC_W, ROAD_LOC_W, PART_ART_D, PART_SEC_D, PART_LOC_D, PART_OG_D, PART_OG_W, PLOT_ART_W, PLOT_SEC_W, PLOT_LOC_W, BLK_ART_NUM_OG_D, BLK_ART_NUM_OG_W, BLK_SEC_NUM_OG_D, BLK_SEC_NUM_OG_W, BLK_LOC_NUM_OG_D, BLK_LOC_NUM_OG_W, PATH_W, OPEN_PERCENT, AMEN_PERCENT, PAVEMENT_W, ADD_TREES, TREE_SPACING, TREE_HEIGHT_START, TREE_HEIGHT_MAX) {

  var __model__ = null;

  /** * **/

  async function exec_cluster1_Zone($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {


    async function exec_cluster1_Zone_node_xlp2w5q5fr_rotSite_($p, pline_art_site_, site_) {

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


    async function exec_cluster1_Zone_node_xlp2w5q5fr_getPerimPlines_($p, site_, road_descr_) {

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


    async function exec_cluster1_Zone_node_xlp2w5q5fr_extendPline_($p, plines_, dist_) {

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

    async function exec_cluster1_Zone_node_xlp2w5q5fr($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
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

      let pline_art_site_ = await exec_cluster1_Zone_node_xlp2w5q5fr_getPerimPlines_($p, site1_, "road_art");
      if ($p.terminated) {
        return mfn.getModel();
      }

      mfn.attrib.Set(null, `pln`, await exec_cluster1_Zone_node_xlp2w5q5fr_rotSite_($p, pline_art_site_, site1_));
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_Zone_node_xlp2w5q5fr_extendPline_($p, pline_art_site_, 200);
      if ($p.terminated) {
        return mfn.getModel();
      }

      mfn.attrib.Set(pline_art_site_, `name`, "art_pline");
    }


    async function exec_cluster1_Zone_node_9ylsi2n679i_getLongestEdge_($p, edge_) {

      let road_ = mfn.attrib.Get(edge_, 'road');

      let edges_same_type_ = [edge_];

      let this_edge_ = edge_;

      let this_road_ = road_;

      while (this_road_ == road_) {

        let force_stop_mobius0_ = 'check_force_stop';

        this_edge_ = await exec_cluster1_Zone_node_9ylsi2n679i_nextEdge_($p, this_edge_);
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

        this_edge_ = await exec_cluster1_Zone_node_9ylsi2n679i_prevEdge_($p, this_edge_);
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


    async function exec_cluster1_Zone_node_9ylsi2n679i_getSiteLongestEdges_($p, site_, road_descr_) {

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


    async function exec_cluster1_Zone_node_9ylsi2n679i_extend_($p, pline_, dist_) {

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


    async function exec_cluster1_Zone_node_9ylsi2n679i_offsetEdge_($p, site_, pline_, dist_) {

      let pline_long_ = await exec_cluster1_Zone_node_9ylsi2n679i_extend_($p, pline_, 1000);
      if ($p.terminated) {
        return mfn.getModel();
      }

      mfn.modify.Offset(pline_long_, -dist_);

      pline_ = await exec_cluster1_Zone_node_9ylsi2n679i_siteTrimPline_($p, site_, pline_long_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      return pline_;
    }


    async function exec_cluster1_Zone_node_9ylsi2n679i_siteTrimPline_($p, site_, pline_) {

      let pline2_ = mfn.poly2d.Boolean(pline_, site_, 'intersect');

      mfn.edit.Delete(pline_, 'delete_selected');

      return pline2_;
    }


    async function exec_cluster1_Zone_node_9ylsi2n679i_getEdge_($p, wire_, posi_) {

      let posis_ = mfn.query.Get('ps', wire_);

      let edges_ = mfn.query.Get('_e', wire_);

      let idx_ = ifn.listFind(posis_, posi_);

      if (idx_ == -1) {

        return null;
      }

      return edges_[pythonList(idx_, edges_.length)];
    }


    async function exec_cluster1_Zone_node_9ylsi2n679i_nextEdge_($p, edge_) {

      let wire_ = mfn.query.Get('_w', edge_);

      let edges_ = mfn.query.Get('_e', wire_);

      let idx_ = ifn.listFind(edges_, edge_);

      let next_idx_ = idx_ + 1;

      if (next_idx_ == ifn.len(edges_)) {

        return edges_[pythonList(0, edges_.length)];
      }

      return edges_[pythonList(next_idx_, edges_.length)];
    }


    async function exec_cluster1_Zone_node_9ylsi2n679i_prevEdge_($p, edge_) {

      let wire_ = mfn.query.Get('_w', edge_);

      let edges_ = mfn.query.Get('_e', wire_);

      let idx_ = ifn.listFind(edges_, edge_);

      let prev_idx_ = idx_ - 1;

      if (prev_idx_ < 0) {

        return edges_[pythonList(-1, edges_.length)];
      }

      return edges_[pythonList(prev_idx_, edges_.length)];
    }


    async function exec_cluster1_Zone_node_9ylsi2n679i_sortAlongX_($p, plines_) {

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


    async function exec_cluster1_Zone_node_9ylsi2n679i_closestEdge_($p, site_, xyz_) {

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

    async function exec_cluster1_Zone_node_9ylsi2n679i($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: left and right plines', '__null__')
      }


      let site1_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', 'main_site')[pythonList(0, 'main_site'.length)];

      let site_art_ = (mfn.query.Filter(mfn.query.Get('_e', null), ['road', null], '==', "road_art"));

      let posis_ = mfn.query.Get('ps', site_art_);

      let edge0_ = await exec_cluster1_Zone_node_9ylsi2n679i_getEdge_($p, site1_, posis_[pythonList(0, posis_.length)]);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let left_edge_ = await exec_cluster1_Zone_node_9ylsi2n679i_prevEdge_($p, edge0_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let left_type_ = mfn.attrib.Get(left_edge_, 'road');

      if (left_type_ != "cold") {

        let left_edge_long_ = await exec_cluster1_Zone_node_9ylsi2n679i_getLongestEdge_($p, left_edge_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        let left_pline_ = await exec_cluster1_Zone_node_9ylsi2n679i_extend_($p, left_edge_long_, 500);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(left_pline_, `name`, "left_pline");

        mfn.attrib.Set(left_pline_, `type`, left_type_);

        mfn.attrib.Set(left_pline_, `site_edge`, left_edge_long_);
      }

      let right_edge_ = await exec_cluster1_Zone_node_9ylsi2n679i_getEdge_($p, site1_, posis_[pythonList(-1, posis_.length)]);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let right_type_ = mfn.attrib.Get(right_edge_, 'road');

      if (right_type_ != "cold") {

        let right_edge_long_ = await exec_cluster1_Zone_node_9ylsi2n679i_getLongestEdge_($p, right_edge_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        let right_pline_ = await exec_cluster1_Zone_node_9ylsi2n679i_extend_($p, right_edge_long_, 500);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(right_pline_, `name`, "right_pline");

        mfn.attrib.Set(right_pline_, `type`, right_type_);

        mfn.attrib.Set(right_pline_, `site_edge`, right_edge_long_);
      }
    }


    async function exec_cluster1_Zone_node_aj9qrn5t0cs_getLongestEdge_($p, edge_) {

      let road_ = mfn.attrib.Get(edge_, 'road');

      let edges_same_type_ = [edge_];

      let this_edge_ = edge_;

      let this_road_ = road_;

      while (this_road_ == road_) {

        this_edge_ = await exec_cluster1_Zone_node_aj9qrn5t0cs_nextEdge_($p, this_edge_);
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

        this_edge_ = await exec_cluster1_Zone_node_aj9qrn5t0cs_prevEdge_($p, this_edge_);
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


    async function exec_cluster1_Zone_node_aj9qrn5t0cs_getSiteLongestEdges_($p, site_, road_descr_) {

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


    async function exec_cluster1_Zone_node_aj9qrn5t0cs_extend_($p, pline_, dist_) {

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


    async function exec_cluster1_Zone_node_aj9qrn5t0cs_offset_($p, sites_, roads_, dist_) {

      let off1_ = mfn.poly2d.OffsetMitre(roads_, dist_, dist_, 'square_end');

      let plines1_ = mfn.make.Polyline(off1_, 'close');

      let plines2_ = mfn.poly2d.Boolean(plines1_, sites_, 'intersect');

      mfn.edit.Delete([off1_, plines1_], 'delete_selected');

      return plines2_;
    }


    async function exec_cluster1_Zone_node_aj9qrn5t0cs_expandEdge_($p, pline_, dist_) {

      let pline_long_ = await exec_cluster1_Zone_node_aj9qrn5t0cs_extend_($p, pline_, 100);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let pgon_ = mfn.poly2d.OffsetMitre(pline_long_, dist_, dist_, 'square_end');

      return pgon_;
    }


    async function exec_cluster1_Zone_node_aj9qrn5t0cs_siteTrimPline_($p, site_, pline_) {

      let pline2_ = mfn.poly2d.Boolean(pline_, site_, 'intersect');

      mfn.edit.Delete(pline_, 'delete_selected');

      return pline2_;
    }


    async function exec_cluster1_Zone_node_aj9qrn5t0cs_getEdge_($p, wire_, posi_) {

      let posis_ = mfn.query.Get('ps', wire_);

      let edges_ = mfn.query.Get('_e', wire_);

      let idx_ = ifn.listFind(posis_, posi_);

      if (idx_ == -1) {

        return null;
      }

      return edges_[pythonList(idx_, edges_.length)];
    }


    async function exec_cluster1_Zone_node_aj9qrn5t0cs_nextEdge_($p, edge_) {

      let wire_ = mfn.query.Get('_w', edge_);

      let edges_ = mfn.query.Get('_e', wire_);

      let idx_ = ifn.listFind(edges_, edge_);

      let next_idx_ = idx_ + 1;

      if (next_idx_ == ifn.len(edges_)) {

        return edges_[pythonList(0, edges_.length)];
      }

      return edges_[pythonList(next_idx_, edges_.length)];
    }


    async function exec_cluster1_Zone_node_aj9qrn5t0cs_prevEdge_($p, edge_) {

      let wire_ = mfn.query.Get('_w', edge_);

      let edges_ = mfn.query.Get('_e', wire_);

      let idx_ = ifn.listFind(edges_, edge_);

      let prev_idx_ = idx_ - 1;

      if (prev_idx_ < 0) {

        return edges_[pythonList(-1, edges_.length)];
      }

      return edges_[pythonList(prev_idx_, edges_.length)];
    }


    async function exec_cluster1_Zone_node_aj9qrn5t0cs_sortAlongX_($p, plines_) {

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


    async function exec_cluster1_Zone_node_aj9qrn5t0cs_closestEdge_($p, edges_, xyz_) {

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


    async function exec_cluster1_Zone_node_aj9qrn5t0cs_transferEdgeAttribs_($p, from_edges_, to_edges_) {

      for (let to_edge_ of to_edges_) {

        let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

        let from_edge_ = await exec_cluster1_Zone_node_aj9qrn5t0cs_closestEdge_($p, from_edges_, cen_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(to_edge_, `road`, mfn.attrib.Get(from_edge_, 'road'));
      }
    }

    async function exec_cluster1_Zone_node_aj9qrn5t0cs($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: zones', '__null__')
      }


      let site1_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', 'main_site')[pythonList(0, 'main_site'.length)];

      let art_pline_ = mfn.query.Filter(mfn.query.Get('pl', null), ['name', null], '==', "art_pline");

      let left_pline_ = mfn.query.Filter(mfn.query.Get('pl', null), ['name', null], '==', "left_pline");

      let right_pline_ = mfn.query.Filter(mfn.query.Get('pl', null), ['name', null], '==', "right_pline");

      let road_art_off_ = await exec_cluster1_Zone_node_aj9qrn5t0cs_offset_($p, site1_, art_pline_, BLK_ART_D_);
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

        road_left_off_ = await exec_cluster1_Zone_node_aj9qrn5t0cs_siteTrimPline_($p, site1_, left_pline_);
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

        road_right_off_ = await exec_cluster1_Zone_node_aj9qrn5t0cs_siteTrimPline_($p, site1_, right_pline_);
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

      await exec_cluster1_Zone_node_aj9qrn5t0cs_transferEdgeAttribs_($p, from_edges_, to_edges_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      mfn.edit.Delete([roads_, zones_], 'keep_selected');
    }


    async function exec_cluster1_Zone_node_z5vxdi5boqk_siteTrimPline_($p, site_, plines_) {

      let new_plines_ = [];

      for (let pline_ of plines_) {

        let pline2_ = mfn.poly2d.Boolean(pline_, site_, 'intersect');

        mfn.edit.Delete(pline_, 'delete_selected');

        mfn.list.Add(new_plines_, pline2_, 'to_end');
      }

      return new_plines_;
    }

    async function exec_cluster1_Zone_node_z5vxdi5boqk($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
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

        x_plines_ = await exec_cluster1_Zone_node_z5vxdi5boqk_siteTrimPline_($p, zone_, x_plines_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(x_plines_, `type`, "road_loc_perp");
      }

      mfn.edit.Delete(x_plines_, 'keep_selected');
    }


    async function exec_cluster1_Zone_node_03lym4l18kqj_siteTrimPline_($p, site_, plines_) {

      let new_plines_ = [];

      for (let pline_ of plines_) {

        let pline2_ = mfn.poly2d.Boolean(pline_, site_, 'intersect');

        mfn.edit.Delete(pline_, 'delete_selected');

        mfn.list.Add(new_plines_, pline2_, 'to_end');
      }

      return new_plines_;
    }


    async function exec_cluster1_Zone_node_03lym4l18kqj_createRoads_($p, bb_all_, zone_, edges_off_) {

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

        new_plines_ = await exec_cluster1_Zone_node_03lym4l18kqj_siteTrimPline_($p, zone_, new_plines_);
        if ($p.terminated) {
          return mfn.getModel();
        }
      }

      return new_plines_;
    }

    async function exec_cluster1_Zone_node_03lym4l18kqj($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
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

        roads_left_ = await exec_cluster1_Zone_node_03lym4l18kqj_createRoads_($p, bb_all_, zone_left_, edges_left_off_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(roads_left_, `type`, "road_loc_para");
      }

      let edges_right_off_ = mfn.query.Filter(mfn.query.Get('_e', zone_right_), ['road', null], '==', "road_loc_right_off");

      let roads_right_ = [];

      if (ifn.len(edges_right_off_) > 0) {

        roads_right_ = await exec_cluster1_Zone_node_03lym4l18kqj_createRoads_($p, bb_all_, zone_right_, edges_right_off_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(roads_right_, `type`, "road_loc_para");
      }

      mfn.edit.Delete([roads_left_, roads_right_], 'keep_selected');
    }


    async function exec_cluster1_Zone_node_vu1s15ynip_siteTrimPline_($p, site_, plines_) {

      let new_plines_ = [];

      for (let pline_ of plines_) {

        let pline2_ = mfn.poly2d.Boolean(pline_, site_, 'intersect');

        mfn.edit.Delete(pline_, 'delete_selected');

        mfn.list.Add(new_plines_, pline2_, 'to_end');
      }

      return new_plines_;
    }

    async function exec_cluster1_Zone_node_vu1s15ynip($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
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

        plines_x_ = await exec_cluster1_Zone_node_vu1s15ynip_siteTrimPline_($p, zone_mid_, plines_x_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(plines_x_, `type`, "road_loc_para_inn");

        plines_y_ = plines_y_.slice(1, -1);

        plines_y_ = await exec_cluster1_Zone_node_vu1s15ynip_siteTrimPline_($p, zone_mid_, plines_y_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(plines_y_, `type`, "road_loc_perp_inn");
      }

      mfn.edit.Delete([plines_x_, plines_y_], 'keep_selected');

      let var_ = mfn.edit.Fuse(mfn.query.Get('pl', null), 0.01);
    }


    async function exec_cluster1_Zone_node_ulnws3299id($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: local ring roads', '__null__')
      }


      mfn.edit.Delete(mfn.query.Get('pl', null), 'keep_selected');
    }


    async function exec_cluster1_Zone_node_93en2qi3yth($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: site', '__null__')
      }


      mfn.edit.Delete(mfn.query.Get('pg', null), 'keep_selected');
    }


    async function exec_cluster1_Zone_node_cpfgl91qjpc_touchingEdge_($p, edges_, xyz_) {

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


    async function exec_cluster1_Zone_node_cpfgl91qjpc_extend_($p, edge_, dist_) {

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


    async function exec_cluster1_Zone_node_cpfgl91qjpc_detachRoadEndsFromColdEdge_($p, site_edges_, roads_) {

      let roads_trimmed_ = [];

      for (let road_ of roads_) {

        let still_exists_ = await exec_cluster1_Zone_node_cpfgl91qjpc_detachFromColdEdge_($p, site_edges_, road_, true);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (still_exists_) {

          still_exists_ = await exec_cluster1_Zone_node_cpfgl91qjpc_detachFromColdEdge_($p, site_edges_, road_, false);
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


    async function exec_cluster1_Zone_node_cpfgl91qjpc_detachFromColdEdge_($p, site_edges_, road_, detach_start_) {

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

      let start_edge_ = await exec_cluster1_Zone_node_cpfgl91qjpc_touchingEdge_($p, site_edges_, xyz_);
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

              await exec_cluster1_Zone_node_cpfgl91qjpc_remDanglingEdges_($p, edges_[pythonList(0, edges_.length)]);
              if ($p.terminated) {
                return mfn.getModel();
              }
            } else {

              await exec_cluster1_Zone_node_cpfgl91qjpc_remDanglingEdges_($p, edges_[pythonList(-1, edges_.length)]);
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


    async function exec_cluster1_Zone_node_cpfgl91qjpc_remDanglingEdges_($p, edge_) {

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


    async function exec_cluster1_Zone_node_cpfgl91qjpc_extendBeyondBoundaryRoad_($p, site_edges_, roads_) {

      for (let road_ of roads_) {

        let posis_ = mfn.query.Get('ps', road_);

        let edges_ = mfn.query.Get('_e', road_);

        let start_edge_ = await exec_cluster1_Zone_node_cpfgl91qjpc_touchingEdge_($p, site_edges_, mfn.attrib.Get(posis_[pythonList(0, posis_.length)], 'xyz'));
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (start_edge_ != null) {

          let start_road_ = mfn.attrib.Get(start_edge_, 'road');

          if (ifn.listFind(["road_art", "road_sec", "road_ter"], start_road_) != -1) {

            await exec_cluster1_Zone_node_cpfgl91qjpc_extend_($p, edges_[pythonList(0, edges_.length)], -20);
            if ($p.terminated) {
              return mfn.getModel();
            }
          }
        }

        let end_edge_ = await exec_cluster1_Zone_node_cpfgl91qjpc_touchingEdge_($p, site_edges_, mfn.attrib.Get(posis_[pythonList(-1, posis_.length)], 'xyz'));
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (end_edge_ != null) {

          let end_road_ = mfn.attrib.Get(end_edge_, 'road');

          if (ifn.listFind(["road_art", "road_sec", "road_ter"], end_road_) != -1) {

            await exec_cluster1_Zone_node_cpfgl91qjpc_extend_($p, edges_[pythonList(-1, edges_.length)], 20);
            if ($p.terminated) {
              return mfn.getModel();
            }
          }
        }
      }
    }


    async function exec_cluster1_Zone_node_cpfgl91qjpc_angDot_($p, vert_) {

      let edges_ = mfn.query.Get('_e', vert_);

      let vecs_ = mfn.calc.Vector(edges_);

      let vec0_ = ifn.vecNorm(vecs_[pythonList(0, vecs_.length)]);

      let vec1_ = ifn.vecNorm(vecs_[pythonList(1, vecs_.length)]);

      return ifn.vecDot(vec0_, vec1_);
    }

    async function exec_cluster1_Zone_node_cpfgl91qjpc($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: roads', '__null__')
      }


      let site1_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', 'main_site')[pythonList(0, 'main_site'.length)];

      let roads_ = mfn.query.Get('pl', null);

      let new_roads_ = mfn.poly2d.Stitch(roads_, 0.1);

      mfn.edit.Delete(roads_, 'delete_selected');

      let site_edges_ = mfn.query.Get('_e', site1_);

      new_roads_ = await exec_cluster1_Zone_node_cpfgl91qjpc_detachRoadEndsFromColdEdge_($p, site_edges_, new_roads_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      new_roads_ = mfn.query.Get('pl', new_roads_);

      await exec_cluster1_Zone_node_cpfgl91qjpc_extendBeyondBoundaryRoad_($p, site_edges_, new_roads_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      new_roads_ = mfn.query.Get('pl', new_roads_);

      let var_ = mfn.edit.Fuse(new_roads_, 5);

      for (let posi_ of mfn.query.Get('ps', new_roads_)) {

        if (ifn.len(mfn.query.Get('_e', posi_)) == 1 && mfn.attrib.Get(posi_, 'type') != "extended") {

          await exec_cluster1_Zone_node_cpfgl91qjpc_remDanglingEdges_($p, mfn.query.Get('_e', posi_));
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

          let dot_ = await exec_cluster1_Zone_node_cpfgl91qjpc_angDot_($p, verts_[pythonList(0, verts_.length)]);
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


    async function exec_cluster1_Zone_node_8joyldrk8dk_transferEdgeAttribs_($p, from_edges_, to_edges_) {

      for (let to_edge_ of to_edges_) {

        let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

        let from_edge_ = await exec_cluster1_Zone_node_8joyldrk8dk_touchingEdge_($p, from_edges_, cen_);
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


    async function exec_cluster1_Zone_node_8joyldrk8dk_touchingEdge_($p, edges_, xyz_) {

      for (let edge_ of edges_) {

        let r_ = mfn.calc.Ray(edge_);

        let d_ = ifn.distance(xyz_, ifn.project(xyz_, r_, 0));

        if (d_ < 1) {

          return edge_;
        }
      }
    }


    async function exec_cluster1_Zone_node_8joyldrk8dk_joinPlines_($p, pline0_, pline1_) {

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


    async function exec_cluster1_Zone_node_8joyldrk8dk__cleanPgonsEdge_($p, pgons_) {

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


    async function exec_cluster1_Zone_node_8joyldrk8dk__cleanPgonsAng_($p, pgons_) {

      for (let pgon_ of pgons_) {

        let del_posis_ = [];

        for (let vert_ of mfn.query.Get('_v', pgon_)) {

          let dot_ = await exec_cluster1_Zone_node_8joyldrk8dk__angDot_($p, vert_);
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


    async function exec_cluster1_Zone_node_8joyldrk8dk__angDot_($p, vert_) {

      let edges_ = mfn.query.Get('_e', vert_);

      let vecs_ = mfn.calc.Vector(edges_);

      let vec0_ = ifn.vecNorm(vecs_[pythonList(0, vecs_.length)]);

      let vec1_ = ifn.vecNorm(vecs_[pythonList(1, vecs_.length)]);

      return ifn.vecDot(vec0_, vec1_);
    }


    async function exec_cluster1_Zone_node_8joyldrk8dk_cleanPgons_($p, pgons_) {

      pgons_ = await exec_cluster1_Zone_node_8joyldrk8dk__cleanPgonsEdge_($p, pgons_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      pgons_ = await exec_cluster1_Zone_node_8joyldrk8dk__cleanPgonsAng_($p, pgons_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      return pgons_;
    }

    async function exec_cluster1_Zone_node_8joyldrk8dk($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
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

            let new_pline_ = await exec_cluster1_Zone_node_8joyldrk8dk_joinPlines_($p, plines_[pythonList(0, plines_.length)], plines_[pythonList(1, plines_.length)]);
            if ($p.terminated) {
              return mfn.getModel();
            }

            mfn.list.Add(bool_roads_, new_pline_, 'to_end');
          }
        }

        bool_roads_ = mfn.query.Get('pl', bool_roads_);

        road_pgons_ = mfn.poly2d.OffsetMitre(mfn.query.Get('pl', bool_roads_), ROAD_LOC_W_ / 2, 100, 'butt_end');

        road_pgons_ = await exec_cluster1_Zone_node_8joyldrk8dk_cleanPgons_($p, road_pgons_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        blocks_ = mfn.poly2d.Boolean(site1_, road_pgons_, 'difference');

        blocks_ = await exec_cluster1_Zone_node_8joyldrk8dk_cleanPgons_($p, blocks_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        road_pgons_ = mfn.poly2d.Boolean(road_pgons_, site1_, 'intersect');

        road_pgons_ = mfn.poly2d.Union(road_pgons_);

        let site_edges_ = mfn.query.Get('_e', site1_);

        for (let block_ of blocks_) {

          await exec_cluster1_Zone_node_8joyldrk8dk_transferEdgeAttribs_($p, site_edges_, mfn.query.Get('_e', block_));
          if ($p.terminated) {
            return mfn.getModel();
          }
        }
      } else {

        blocks_ = mfn.make.Copy(site1_, null);

        let site_edges_ = mfn.query.Get('_e', site1_);

        await exec_cluster1_Zone_node_8joyldrk8dk_transferEdgeAttribs_($p, site_edges_, mfn.query.Get('_e', blocks_));
        if ($p.terminated) {
          return mfn.getModel();
        }
      }

      mfn.edit.Delete([roads_, blocks_, road_pgons_], 'keep_selected');

      mfn.attrib.Set(blocks_, `type`, "block");

      mfn.attrib.Set(road_pgons_, `type`, "road_loc");
    }


    async function exec_cluster1_Zone_node_8tqxj3gvnfc($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: undo rot', '__null__')
      }


      let ents_ = [mfn.query.Get('pl', null), mfn.query.Get('pg', null)];

      mfn.modify.XForm(ents_, JSON.parse(JSON.stringify(ifn.XY)), mfn.attrib.Get(null, 'pln'));
    }


    async function exec_cluster1_Zone_node_9bhi1flny5t($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: End', '__null__')
      }


      return null;
    }

    var merged;
    let ssid_exec_cluster1_Zone_node_8ylatkqudwy = mfn.model.snapshotGetActive();

    let ssid_exec_cluster1_Zone_node_xlp2w5q5fr = mfn.model.snapshotNext([ssid_exec_cluster1_Zone_node_8ylatkqudwy]);

    await exec_cluster1_Zone_node_xlp2w5q5fr($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_cluster1_Zone_node_9ylsi2n679i = mfn.model.snapshotNext([ssid_exec_cluster1_Zone_node_xlp2w5q5fr]);

    await exec_cluster1_Zone_node_9ylsi2n679i($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_cluster1_Zone_node_aj9qrn5t0cs = ssid_exec_cluster1_Zone_node_9ylsi2n679i;

    await exec_cluster1_Zone_node_aj9qrn5t0cs($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_cluster1_Zone_node_z5vxdi5boqk = mfn.model.snapshotNext([ssid_exec_cluster1_Zone_node_aj9qrn5t0cs]);

    await exec_cluster1_Zone_node_z5vxdi5boqk($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_cluster1_Zone_node_03lym4l18kqj = mfn.model.snapshotNext([ssid_exec_cluster1_Zone_node_aj9qrn5t0cs]);

    await exec_cluster1_Zone_node_03lym4l18kqj($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_cluster1_Zone_node_vu1s15ynip = mfn.model.snapshotNext([ssid_exec_cluster1_Zone_node_aj9qrn5t0cs]);

    await exec_cluster1_Zone_node_vu1s15ynip($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_cluster1_Zone_node_ulnws3299id = mfn.model.snapshotNext([ssid_exec_cluster1_Zone_node_aj9qrn5t0cs]);

    await exec_cluster1_Zone_node_ulnws3299id($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_cluster1_Zone_node_93en2qi3yth = mfn.model.snapshotNext([ssid_exec_cluster1_Zone_node_xlp2w5q5fr]);

    await exec_cluster1_Zone_node_93en2qi3yth($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_cluster1_Zone_node_cpfgl91qjpc = mfn.model.snapshotNext([ssid_exec_cluster1_Zone_node_z5vxdi5boqk, ssid_exec_cluster1_Zone_node_03lym4l18kqj, ssid_exec_cluster1_Zone_node_vu1s15ynip, ssid_exec_cluster1_Zone_node_ulnws3299id, ssid_exec_cluster1_Zone_node_93en2qi3yth]);

    await exec_cluster1_Zone_node_cpfgl91qjpc($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_cluster1_Zone_node_8joyldrk8dk = ssid_exec_cluster1_Zone_node_cpfgl91qjpc;

    await exec_cluster1_Zone_node_8joyldrk8dk($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_cluster1_Zone_node_8tqxj3gvnfc = ssid_exec_cluster1_Zone_node_8joyldrk8dk;

    await exec_cluster1_Zone_node_8tqxj3gvnfc($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);

    let ssid_exec_cluster1_Zone_node_9bhi1flny5t = mfn.model.snapshotNext([ssid_exec_cluster1_Zone_node_8ylatkqudwy, ssid_exec_cluster1_Zone_node_8tqxj3gvnfc]);

    return await exec_cluster1_Zone_node_9bhi1flny5t($p, SITE_, BLK_ART_D_, BLK_ART_W_, BLK_SEC_D_, BLK_SEC_W_, BLK_LOC_D_, BLK_LOC_W_, ROAD_LOC_W_);
  }

  async function exec_cluster1($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {

    async function exec_cluster1_node_6w5xaj6iu9j($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: import', '__null__')
      }


      let coll_city_ = await mfn.io.Import(IN_MODEL_, 'sim');

      mfn.attrib.Add('pg', 'string', 'cluster_id');
    }


    async function exec_cluster1_node_f9mnqqhlrz($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: warm blocks', '__null__')
      }


      let all_blocks_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', "block");

      let warm_blocks_ = mfn.query.Filter(all_blocks_, ['block_type', null], '!=', "cold");

      mfn.edit.Delete(mfn.query.Get('pg', warm_blocks_), 'keep_selected');
    }


    async function exec_cluster1_node_8u13jv4q514_getCentralBlock_($p, blocks_, site_cen_) {

      let cen_block_ = null;

      let min_dist_ = Infinity;

      for (let block_ of blocks_) {

        let block_cen_ = mfn.calc.Centroid(block_, 'center_of_mass');

        let dist_ = ifn.distance(site_cen_, block_cen_);

        if (dist_ < min_dist_) {

          min_dist_ = dist_;

          cen_block_ = block_;
        }
      }

      return cen_block_;
    }

    async function exec_cluster1_node_8u13jv4q514($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: central block', '__null__')
      }


      let warm_blocks_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', "block");

      let warm_not_art_blocks_ = mfn.query.Filter(warm_blocks_, ['block_type', null], '!=', "art");

      for (let site_name_ of mfn.attrib.Get(null, 'site_names')) {

        let site_cen_ = mfn.attrib.Get(null, site_name_ + "_cen");

        let site_area_ = mfn.attrib.Get(null, site_name_ + "_area");

        let open_req_area_ = site_area_ * (OPEN_PERCENT_ / 100);

        let amen_req_area_ = site_area_ * (AMEN_PERCENT_ / 100);

        let req_area_ = open_req_area_ + amen_req_area_;

        let blocks_ = mfn.query.Filter(warm_not_art_blocks_, ['site', null], '==', site_name_);

        let cen_block_ = await exec_cluster1_node_8u13jv4q514_getCentralBlock_($p, blocks_, site_cen_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (cen_block_ == null) {

          continue;
        }

        let cen_block_area_ = mfn.calc.Area(cen_block_);

        mfn.visualize.Color(cen_block_, [0, 1, 0]);

        mfn.attrib.Set(null, site_name_ + "_cen_block_id", mfn.attrib.Get(cen_block_, 'block_id'), 'one_value');
      }
    }


    async function exec_cluster1_node_v300sbtwm6p_createExtendedRoads_($p, block_) {

      let roads_art_ = await exec_cluster1_node_v300sbtwm6p__getPerimPlines_($p, block_, "road_art");
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_node_v300sbtwm6p__extendPline_($p, roads_art_, 100);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let roads_sec_ = await exec_cluster1_node_v300sbtwm6p__getPerimPlines_($p, block_, "road_sec");
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_node_v300sbtwm6p__extendPline_($p, roads_sec_, 100);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let roads_loc_ = await exec_cluster1_node_v300sbtwm6p__getPerimPlines_($p, block_, "road_loc");
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_node_v300sbtwm6p__extendPline_($p, roads_loc_, 100);
      if ($p.terminated) {
        return mfn.getModel();
      }

      return [roads_art_, roads_sec_, roads_loc_];
    }


    async function exec_cluster1_node_v300sbtwm6p_createOffGridArea_($p, block_, roads_) {

      let roads_art_ = roads_[pythonList(0, roads_.length)];

      let roads_sec_ = roads_[pythonList(1, roads_.length)];

      let roads_loc_ = roads_[pythonList(2, roads_.length)];

      let off_art_ = mfn.poly2d.OffsetMitre(roads_art_, PART_ART_D_, 100, 'butt_end');

      let off_sec_ = mfn.poly2d.OffsetMitre(roads_sec_, PART_SEC_D_, 100, 'butt_end');

      let off_loc_ = mfn.poly2d.OffsetMitre(roads_loc_, PART_LOC_D_, 100, 'butt_end');

      let off_grids_ = mfn.poly2d.Boolean(block_, [off_art_, off_sec_, off_loc_], 'difference');

      mfn.edit.Delete([off_art_, off_sec_, off_loc_], 'delete_selected');

      off_grids_ = await exec_cluster1_node_v300sbtwm6p__cleanPgonsEdge_($p, off_grids_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      off_grids_ = await exec_cluster1_node_v300sbtwm6p__cleanPgonsAng_($p, off_grids_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      if (ifn.len(off_grids_) == 0) {

        return null;
      }

      let off_grid_ = off_grids_[pythonList(0, off_grids_.length)];

      if (ifn.len(off_grids_) > 1) {

        let areas_ = mfn.calc.Area(off_grids_);

        let sorted_ = ifn.listSort(off_grids_, areas_);

        off_grid_ = sorted_[pythonList(-1, sorted_.length)];

        for (let small_off_grid_ of off_grids_) {

          if (small_off_grid_ != off_grid_) {

            mfn.edit.Delete(small_off_grid_, 'delete_selected');
          }
        }
      }

      mfn.visualize.Color(off_grid_, [0.7, 0.7, 1]);

      mfn.attrib.Set(off_grid_, `block_type`, mfn.attrib.Get(block_, 'block_type'));

      return off_grid_;
    }


    async function exec_cluster1_node_v300sbtwm6p_createBlockPartsFromOffGrid_($p, block_, off_grid_) {

      let outer_edges_ = mfn.query.Get('_e', block_);

      let outer_rays_ = mfn.calc.Ray(outer_edges_);

      let corners1_ = [];

      let side_posis_ = [];

      let perv_vert_ = null;

      for (let vert_ of mfn.query.Get('_v', off_grid_)) {

        let ang_ = await exec_cluster1_node_v300sbtwm6p__vertAng_($p, vert_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (ang_ < 155) {

          let edges_ = mfn.query.Get('_e', vert_);

          let er0_ = mfn.calc.Ray(edges_[pythonList(0, edges_.length)]);

          let er1_ = mfn.calc.Ray(edges_[pythonList(1, edges_.length)]);

          let v0_perp_ = ifn.vecSetLen([er0_[pythonList(1, er0_.length)][pythonList(1, er0_[pythonList(1, er0_.length)].length)], -er0_[pythonList(1, er0_.length)][pythonList(0, er0_[pythonList(1, er0_.length)].length)], 0], 1000);

          let v1_perp_ = ifn.vecSetLen([er1_[pythonList(1, er1_.length)][pythonList(1, er1_[pythonList(1, er1_.length)].length)], -er1_[pythonList(1, er1_.length)][pythonList(0, er1_[pythonList(1, er1_.length)].length)], 0], 1000);

          let vm_perp_ = ifn.vecSetLen(ifn.vecSum(v0_perp_, v1_perp_), 1000);

          let r0_perp_ = ifn.rayMake(er1_[pythonList(0, er1_.length)], v0_perp_);

          let r1_perp_ = ifn.rayMake(er1_[pythonList(0, er1_.length)], v1_perp_);

          let rm_perp_ = ifn.rayMake(er1_[pythonList(0, er1_.length)], ifn.vecSum(v0_perp_, v1_perp_));

          let ori_posi_ = mfn.query.Get('ps', vert_)[pythonList(0, mfn.query.Get('ps', vert_).length)];

          let xyz0_ = ifn.vecAdd(er1_[pythonList(0, er1_.length)], v0_perp_);

          let xyz1_ = ifn.vecAdd(er1_[pythonList(0, er1_.length)], v1_perp_);

          let xyzm_ = ifn.vecAdd(er1_[pythonList(0, er1_.length)], vm_perp_);

          let posis_ = mfn.make.Position([xyz0_, xyzm_, xyz1_]);

          let corner_ = mfn.make.Polygon(ifn.listJoin(posis_, ori_posi_));

          mfn.list.Add(corners1_, corner_, 'to_end');

          mfn.list.Add(side_posis_, [ori_posi_, posis_[pythonList(0, posis_.length)], edges_], 'to_end');

          mfn.list.Add(side_posis_, [ori_posi_, posis_[pythonList(-1, posis_.length)], edges_], 'to_end');
        }
      }

      let sides1_ = [];

      let side_posis_len_ = ifn.len(side_posis_);

      for (let i_ of ifn.range(1, side_posis_len_, 2)) {

        let a_ = side_posis_[pythonList(i_ % side_posis_len_, side_posis_.length)];

        let b_ = side_posis_[pythonList((i_ + 1) % side_posis_len_, side_posis_.length)];

        let perim_posis_ = await exec_cluster1_node_v300sbtwm6p__getPosisBetweenEdges_($p, a_[pythonList(2, a_.length)][pythonList(1, a_[pythonList(2, a_.length)].length)], b_[pythonList(2, b_.length)][pythonList(0, b_[pythonList(2, b_.length)].length)], true);
        if ($p.terminated) {
          return mfn.getModel();
        }

        let posis_side_ = ifn.listFlat([a_[pythonList(0, a_.length)], a_[pythonList(1, a_.length)], b_[pythonList(1, b_.length)], b_[pythonList(0, b_.length)], ifn.listRev(perim_posis_)]);

        let side_ = mfn.make.Polygon(posis_side_);

        mfn.list.Add(sides1_, side_, 'to_end');
      }

      let corner_parts_ = mfn.poly2d.Boolean(corners1_, block_, 'intersect');

      let side_parts_ = mfn.poly2d.Boolean(sides1_, block_, 'intersect');

      mfn.visualize.Color(corner_parts_, [1, 0.7, 0.7]);

      mfn.visualize.Color(side_parts_, [0.7, 1, 0.7]);

      mfn.attrib.Set(side_parts_, `type`, "side");

      let side_parts_e_ = mfn.query.Get('_e', side_parts_);

      let corner_parts_e_ = mfn.query.Get('_e', corner_parts_);

      let off_grid_e_ = mfn.query.Get('_e', off_grid_);

      await exec_cluster1_node_v300sbtwm6p__setAttribs_($p, block_, ifn.listFlat([corner_parts_, side_parts_]), off_grid_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_node_v300sbtwm6p__transferEdgeAttribsTouchingPart_($p, corner_parts_e_, side_parts_e_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_node_v300sbtwm6p__transferEdgeAttribsTouchingPart_($p, side_parts_e_, corner_parts_e_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_node_v300sbtwm6p__transferEdgeAttribsTouchingPart_($p, off_grid_e_, side_parts_e_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_node_v300sbtwm6p__transferEdgeAttribsTouchingPart_($p, side_parts_e_, off_grid_e_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      mfn.edit.Delete([corners1_, sides1_, block_], 'delete_selected');

      return [corner_parts_, side_parts_, off_grid_];
    }


    async function exec_cluster1_node_v300sbtwm6p__isectRays_($p, ray_, rays_) {

      for (let i_ of ifn.range(ifn.len(rays_))) {

        let isect_ = ifn.intersect(ray_, rays_[pythonList(i_, rays_.length)], 0);

        if (isect_ != null) {

          return [isect_, i_];
        }
      }

      return null;
    }


    async function exec_cluster1_node_v300sbtwm6p__cleanPgonsEdge_($p, pgons_) {

      for (let pgon_ of pgons_) {

        let del_posis_ = [];

        for (let edge_ of mfn.query.Get('_e', pgon_)) {

          let length_ = mfn.calc.Length(edge_);

          if (length_ < 0.1) {

            let posis_ = mfn.query.Get('ps', edge_);

            mfn.list.Add(del_posis_, posis_[pythonList(0, posis_.length)], 'to_end');
          }
        }
      }

      return mfn.query.Get('pg', pgons_);
    }


    async function exec_cluster1_node_v300sbtwm6p__cleanPgonsAng_($p, pgons_) {

      for (let pgon_ of pgons_) {

        let del_posis_ = [];

        for (let vert_ of mfn.query.Get('_v', pgon_)) {

          let dot_ = await exec_cluster1_node_v300sbtwm6p__angDot_($p, vert_);
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


    async function exec_cluster1_node_v300sbtwm6p__angDot_($p, vert_) {

      let edges_ = mfn.query.Get('_e', vert_);

      let vecs_ = mfn.calc.Vector(edges_);

      let vec0_ = ifn.vecNorm(vecs_[pythonList(0, vecs_.length)]);

      let vec1_ = ifn.vecNorm(vecs_[pythonList(1, vecs_.length)]);

      return ifn.vecDot(vec0_, vec1_);
    }


    async function exec_cluster1_node_v300sbtwm6p__vertAng_($p, vert_) {

      let edges_ = mfn.query.Get('_e', vert_);

      if (ifn.len(edges_) == 1) {

        return 0;
      }

      let vecs_ = mfn.calc.Vector(edges_);

      let vec0_ = ifn.vecRev(ifn.vecNorm(vecs_[pythonList(0, vecs_.length)]));

      let vec1_ = ifn.vecNorm(vecs_[pythonList(1, vecs_.length)]);

      return ifn.radToDeg(ifn.vecAng2(vec1_, vec0_, [0, 0, 1]));
    }


    async function exec_cluster1_node_v300sbtwm6p__getPerimPlines_($p, site_, road_descr_) {

      let posis_ = [];

      for (let edge_ of mfn.query.Get('_e', site_)) {

        if (mfn.attrib.Get(edge_, 'road') == road_descr_) {

          let start_end_ = mfn.query.Get('ps', edge_);

          if (ifn.len(posis_) == 0 || start_end_[pythonList(0, start_end_.length)] != posis_[pythonList(-1, posis_.length)][pythonList(-1, posis_[pythonList(-1, posis_.length)].length)]) {

            mfn.list.Add(posis_, start_end_, 'to_end');
          } else {

            mfn.list.Add(posis_[pythonList(-1, posis_.length)], start_end_[pythonList(1, start_end_.length)], 'to_end');
          }
        }
      }

      if (ifn.len(posis_) == 0) {

        return [];
      }

      if (ifn.len(posis_) > 1 && posis_[pythonList(-1, posis_.length)][pythonList(-1, posis_[pythonList(-1, posis_.length)].length)] == posis_[pythonList(0, posis_.length)][pythonList(0, posis_[pythonList(0, posis_.length)].length)]) {

        let first_list_ = ifn.listJoin(posis_[pythonList(-1, posis_.length)], posis_[pythonList(0, posis_.length)].slice(1));

        posis_[pythonList(0, posis_.length)] = first_list_;

        posis_ = posis_.slice(0, -1);
      }

      let site_plines_ = mfn.make.Polyline(posis_, 'open');

      return site_plines_;
    }


    async function exec_cluster1_node_v300sbtwm6p__setAttribs_($p, block_, parts_, off_grid_) {

      for (let part_ of parts_) {

        await exec_cluster1_node_v300sbtwm6p__transferEdgeAttribsTouching_($p, mfn.query.Get('_e', block_), mfn.query.Get('_e', part_));
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (mfn.attrib.Get(part_, 'type') == "side") {

          await exec_cluster1_node_v300sbtwm6p__addAttribs_($p, part_, true);
          if ($p.terminated) {
            return mfn.getModel();
          }
        } else {

          await exec_cluster1_node_v300sbtwm6p__addAttribs_($p, part_, false);
          if ($p.terminated) {
            return mfn.getModel();
          }
        }
      }

      if (off_grid_ == null) {

        let var_ = mfn.edit.Fuse(parts_, 0.001);
      } else {

        await exec_cluster1_node_v300sbtwm6p__addAttribs_($p, off_grid_, false);
        if ($p.terminated) {
          return mfn.getModel();
        }

        let var_ = mfn.edit.Fuse(ifn.listFlat([parts_, off_grid_]), 0.001);
      }
    }


    async function exec_cluster1_node_v300sbtwm6p__touchingEdge_($p, edges_, xyz_) {

      for (let edge_ of edges_) {

        let r_ = mfn.calc.Ray(edge_);

        let d_ = ifn.distance(xyz_, ifn.project(xyz_, r_, 0));

        if (d_ < 1) {

          return edge_;
        }
      }

      return null;
    }


    async function exec_cluster1_node_v300sbtwm6p__extendPline_($p, plines_, dist_) {

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


    async function exec_cluster1_node_v300sbtwm6p__addAttribs_($p, part_, is_side_part_) {

      let road_types_ = [];

      for (let edge_ of mfn.query.Get('_e', part_)) {

        let road_type_ = mfn.attrib.Get(edge_, 'road');

        if (road_type_ != undefined) {

          mfn.list.Add(road_types_, road_type_, 'to_end');
        }
      }

      let check_ = road_types_;

      mfn.attrib.Set(part_, `class`, "part");

      let block_types_dict_ = {
        "road_art": "art",
        "road_sec": "sec",
        "road_ter": "ter",
        "road_loc": "loc"
      };

      let cats_dict_ = {
        "road_art": 1,
        "road_sec": 2,
        "road_ter": 3,
        "road_loc": 4
      };

      if (ifn.len(road_types_) == 0) {

        mfn.attrib.Set(part_, `type`, "off_grid");
      } else {
        if (ifn.len(road_types_) == 1) {

          let road_type_ = road_types_[pythonList(0, road_types_.length)];

          let block_type_ = ifn.string(block_types_dict_[pythonList(road_type_, block_types_dict_.length)]);

          mfn.attrib.Set(part_, `type`, block_type_);
        } else {

          let data_ = [];

          let data_sort_ = [];

          for (let road_type_ of road_types_) {

            mfn.list.Add(data_, [road_type_, cats_dict_[pythonList(road_type_, cats_dict_.length)]], 'to_end');

            mfn.list.Add(data_sort_, cats_dict_[pythonList(road_type_, cats_dict_.length)], 'to_end');
          }

          let sorted_data_ = ifn.listSort(data_, data_sort_);

          let road_type0_ = sorted_data_[pythonList(0, sorted_data_.length)][pythonList(0, sorted_data_[pythonList(0, sorted_data_.length)].length)];

          let block_type0_ = ifn.string(block_types_dict_[pythonList(road_type0_, block_types_dict_.length)]);

          let road_type1_ = sorted_data_[pythonList(1, sorted_data_.length)][pythonList(0, sorted_data_[pythonList(1, sorted_data_.length)].length)];

          let block_type1_ = ifn.string(block_types_dict_[pythonList(road_type1_, block_types_dict_.length)]);

          if (is_side_part_) {

            mfn.attrib.Set(part_, `type`, block_type0_);
          } else {

            mfn.attrib.Set(part_, `type`, block_type0_ + "_" + block_type1_);
          }
        }
      }
    }


    async function exec_cluster1_node_v300sbtwm6p__transferEdgeAttribsTouching_($p, from_edges_, to_edges_) {

      for (let to_edge_ of to_edges_) {

        let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

        let from_edge_ = await exec_cluster1_node_v300sbtwm6p__touchingEdge_($p, from_edges_, cen_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (from_edge_ != null) {

          mfn.attrib.Set(to_edge_, `road`, mfn.attrib.Get(from_edge_, 'road'));
        }
      }
    }


    async function exec_cluster1_node_v300sbtwm6p__transferEdgeAttribsTouchingPart_($p, from_edges_, to_edges_) {

      for (let to_edge_ of to_edges_) {

        if (mfn.attrib.Get(to_edge_, 'road') == undefined) {

          let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

          let from_edge_ = await exec_cluster1_node_v300sbtwm6p__touchingEdge_($p, from_edges_, cen_);
          if ($p.terminated) {
            return mfn.getModel();
          }

          if (from_edge_ != null) {

            mfn.attrib.Set(to_edge_, `road`, mfn.attrib.Get(mfn.query.Get('pg', from_edge_)[pythonList(0, mfn.query.Get('pg', from_edge_).length)], 'type'));
          }
        }
      }
    }


    async function exec_cluster1_node_v300sbtwm6p__getPosisBetweenEdges_($p, edge0_, edge1_, forwards_) {

      let idx_ = 1;

      if (!forwards_) {

        idx_ = 0;
      }

      let posis_ = [];

      let next_edge_ = edge0_;

      while (next_edge_ != edge1_) {

        let edge_posis_ = mfn.query.Get('ps', next_edge_);

        mfn.list.Add(posis_, edge_posis_[pythonList(1, edge_posis_.length)], 'to_end');

        next_edge_ = mfn.query.Get('_e', mfn.query.Get('_v', next_edge_)[pythonList(idx_, mfn.query.Get('_v', next_edge_).length)])[pythonList(idx_, mfn.query.Get('_e', mfn.query.Get('_v', next_edge_)[pythonList(idx_, mfn.query.Get('_v', next_edge_).length)]).length)];
      }

      return posis_;
    }

    async function exec_cluster1_node_v300sbtwm6p($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: ', '__null__')
      }


      let blocks_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', "block");

      for (let block_ of blocks_) {

        let roads_ = await exec_cluster1_node_v300sbtwm6p_createExtendedRoads_($p, block_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        let off_grid_ = await exec_cluster1_node_v300sbtwm6p_createOffGridArea_($p, block_, roads_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (off_grid_ != null) {

          let area_ = mfn.calc.Area(off_grid_);

          if (area_ < (PART_OG_D_ * PART_OG_W_ * 0.5)) {

            mfn.edit.Delete(off_grid_, 'delete_selected');

            mfn.attrib.Set(block_, `class`, "block");
          } else {

            let site_ = mfn.attrib.Get(block_, 'site');

            let block_id_ = mfn.attrib.Get(block_, 'block_id');

            let block_type_ = mfn.attrib.Get(block_, 'block_type');

            let parts_ = await exec_cluster1_node_v300sbtwm6p_createBlockPartsFromOffGrid_($p, block_, off_grid_);
            if ($p.terminated) {
              return mfn.getModel();
            }

            mfn.attrib.Set(parts_, `class`, "part");

            mfn.attrib.Set(parts_, `site`, site_);

            mfn.attrib.Set(parts_, `block_type`, block_type_);

            mfn.attrib.Set(parts_, `block_id`, block_id_);
          }
        }

        mfn.edit.Delete(roads_, 'delete_selected');
      }
    }


    async function exec_cluster1_node_9bpb8dz1lre_createParts_($p, block_, depths_dict_) {

      let plines_art_ = await exec_cluster1_node_9bpb8dz1lre__getPerimPlines_($p, block_, "road_art");
      if ($p.terminated) {
        return mfn.getModel();
      }

      let plines_sec_ = await exec_cluster1_node_9bpb8dz1lre__getPerimPlines_($p, block_, "road_sec");
      if ($p.terminated) {
        return mfn.getModel();
      }

      if (ifn.len(plines_art_) > 0 && ifn.len(plines_sec_) == 0) {

        let art_off_ = mfn.poly2d.OffsetMitre(plines_art_, PART_ART_D_, PART_ART_D_, 'square_end');

        let art_bool1_ = mfn.poly2d.Boolean(art_off_, block_, 'intersect');

        let art_bool2_ = mfn.poly2d.Boolean(block_, art_off_, 'difference');

        let result_ = await exec_cluster1_node_9bpb8dz1lre__trimLoc_($p, block_, [art_bool1_, art_bool2_]);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.edit.Delete([art_off_, plines_art_], 'delete_selected');

        return result_;
      } else {
        if (ifn.len(plines_art_) == 0 && ifn.len(plines_sec_) > 0) {

          let sec_off_ = mfn.poly2d.OffsetMitre(plines_sec_, PART_SEC_D_, PART_SEC_D_, 'square_end');

          let sec_bool1_ = mfn.poly2d.Boolean(sec_off_, block_, 'intersect');

          let sec_bool2_ = mfn.poly2d.Boolean(block_, sec_off_, 'difference');

          let result_ = await exec_cluster1_node_9bpb8dz1lre__trimLoc_($p, block_, [sec_bool1_, sec_bool2_]);
          if ($p.terminated) {
            return mfn.getModel();
          }

          mfn.edit.Delete([sec_off_, plines_sec_], 'delete_selected');

          return result_;
        } else {
          if (ifn.len(plines_art_) > 0 && ifn.len(plines_sec_) > 0) {

            let art_off_ = mfn.poly2d.OffsetMitre(plines_art_, PART_ART_D_, PART_ART_D_, 'square_end');

            let art_bool1_ = mfn.poly2d.Boolean(art_off_, block_, 'intersect');

            let art_bool2_ = mfn.poly2d.Boolean(block_, art_off_, 'difference');

            let sec_off_ = mfn.poly2d.OffsetMitre(plines_sec_, PART_SEC_D_, PART_SEC_D_, 'square_end');

            let sec_bool1_ = mfn.poly2d.Boolean(sec_off_, block_, 'intersect');

            let sec_bool2_ = mfn.poly2d.Boolean(block_, sec_off_, 'difference');

            let art_corners_ = mfn.poly2d.Boolean(art_bool1_, sec_off_, 'intersect');

            let other_corners_ = mfn.poly2d.Boolean(art_bool2_, sec_off_, 'intersect');

            let art_bool1_trim_ = mfn.poly2d.Boolean(art_bool1_, sec_off_, 'difference');

            let art_bool2_trim_ = mfn.poly2d.Boolean(art_bool2_, sec_off_, 'difference');

            mfn.edit.Delete([art_off_, plines_art_, sec_off_, plines_sec_, art_bool1_, art_bool2_, sec_bool1_, sec_bool2_], 'delete_selected');

            mfn.visualize.Color(art_bool2_trim_, [1, 0, 0]);

            let art_bool2_trim_areas_ = mfn.calc.Area(art_bool2_trim_);

            if (ifn.sum(art_bool2_trim_areas_) > (PART_LOC_D_ * PART_LOC_D_ * 4)) {

              let result_ = await exec_cluster1_node_9bpb8dz1lre__trimLoc_($p, block_, [art_corners_, other_corners_, art_bool1_trim_, art_bool2_trim_]);
              if ($p.terminated) {
                return mfn.getModel();
              }

              return result_;
            } else {

              return ifn.listFlat([art_corners_, other_corners_, art_bool1_trim_, art_bool2_trim_]);
            }
          }
        }
      }

      return [];
    }


    async function exec_cluster1_node_9bpb8dz1lre_transferEdgeAttribsTouching_($p, from_edges_, to_edges_) {

      for (let to_edge_ of to_edges_) {

        let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

        let from_edge_ = await exec_cluster1_node_9bpb8dz1lre__touchingEdge_($p, from_edges_, cen_);
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


    async function exec_cluster1_node_9bpb8dz1lre_addTypeAttribsCrvRoads_($p, part_) {

      let part_edges_ = mfn.query.Get('_e', part_);

      let road_types_ = [];

      let art_ = await exec_cluster1_node_9bpb8dz1lre__getCrvRoadEdgeTypes_($p, part_edges_, "road_art");
      if ($p.terminated) {
        return mfn.getModel();
      }

      let sec_ = await exec_cluster1_node_9bpb8dz1lre__getCrvRoadEdgeTypes_($p, part_edges_, "road_sec");
      if ($p.terminated) {
        return mfn.getModel();
      }

      let loc_ = await exec_cluster1_node_9bpb8dz1lre__getCrvRoadEdgeTypes_($p, part_edges_, "road_loc");
      if ($p.terminated) {
        return mfn.getModel();
      }

      mfn.list.Add(road_types_, ifn.listFlat([art_, sec_, loc_]), 'extend_end');

      mfn.attrib.Set(part_, `class`, "part");

      let block_types_dict_ = {
        "road_art": "art",
        "road_sec": "sec",
        "road_ter": "ter",
        "road_loc": "loc"
      };

      let cats_dict_ = {
        "road_art": 1,
        "road_sec": 2,
        "road_ter": 3,
        "road_loc": 4
      };

      if (ifn.len(road_types_) == 0) {

        mfn.attrib.Set(part_, `type`, "off_grid");
      } else {
        if (ifn.len(road_types_) == 1) {

          let road_type_ = road_types_[pythonList(0, road_types_.length)];

          let block_type_ = ifn.string(block_types_dict_[pythonList(road_type_, block_types_dict_.length)]);

          mfn.attrib.Set(part_, `type`, block_type_);
        } else {

          let road_type0_ = road_types_[pythonList(0, road_types_.length)];

          let block_type0_ = ifn.string(block_types_dict_[pythonList(road_type0_, block_types_dict_.length)]);

          let road_type1_ = road_types_[pythonList(1, road_types_.length)];

          let block_type1_ = ifn.string(block_types_dict_[pythonList(road_type1_, block_types_dict_.length)]);

          mfn.attrib.Set(part_, `type`, block_type0_ + "_" + block_type1_);
        }
      }
    }


    async function exec_cluster1_node_9bpb8dz1lre_transferEdgeAttribsBtwTouchingParts_($p, parts_) {

      let edges_ = mfn.query.Get('_e', parts_);

      for (let to_edge_ of edges_) {

        if (mfn.attrib.Get(to_edge_, 'road') == undefined) {

          let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

          let idx_ = ifn.listFind(edges_, to_edge_);

          let from_edges_ = ifn.listJoin(edges_.slice(0, idx_), edges_.slice(idx_ + 1));

          let from_edge_ = await exec_cluster1_node_9bpb8dz1lre__touchingEdge_($p, from_edges_, cen_);
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


    async function exec_cluster1_node_9bpb8dz1lre_copyAttribs_($p, pgon_from_, pgons_to_) {

      let block_id_ = mfn.attrib.Get(pgon_from_, 'block_id');

      let block_type_ = mfn.attrib.Get(pgon_from_, 'block_type');

      let site_ = mfn.attrib.Get(pgon_from_, 'site');

      for (let pgon_to_ of pgons_to_) {

        mfn.attrib.Set(pgon_to_, `block_id`, block_id_);

        mfn.attrib.Set(pgon_to_, `block_type`, block_type_);

        mfn.attrib.Set(pgon_to_, `site`, site_);
      }
    }


    async function exec_cluster1_node_9bpb8dz1lre__getCrvRoadEdgeTypes_($p, part_edges_, road_descr_) {

      let edges_ = [[]];

      for (let edge_ of part_edges_) {

        let edge_road_ = mfn.attrib.Get(edge_, 'road');

        if (edge_road_ == road_descr_) {

          mfn.list.Add(edges_[pythonList(-1, edges_.length)], edge_, 'to_end');
        } else {

          if (ifn.len(edges_[pythonList(-1, edges_.length)]) != 0) {

            mfn.list.Add(edges_, [], 'to_end');
          }
        }
      }

      if (ifn.len(edges_) == 2) {

        edges_ = ifn.listJoin(edges_[pythonList(1, edges_.length)], edges_[pythonList(0, edges_.length)]);
      } else {

        edges_ = edges_[pythonList(0, edges_.length)];
      }

      if (ifn.len(edges_) == 0) {

        return [];
      }

      if (ifn.len(edges_) == 1) {

        return [road_descr_];
      }

      let road_types_ = [road_descr_];

      for (let i_ of ifn.range(1, ifn.len(edges_))) {

        let ang_ = await exec_cluster1_node_9bpb8dz1lre__vertAng_($p, mfn.query.Get('_v', edges_[pythonList(i_, edges_.length)])[pythonList(0, mfn.query.Get('_v', edges_[pythonList(i_, edges_.length)]).length)]);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (ifn.abs(ang_ - 180) > 45) {

          mfn.list.Add(road_types_, road_descr_, 'to_end');
        }
      }

      return road_types_;
    }


    async function exec_cluster1_node_9bpb8dz1lre__vertAng_($p, vert_) {

      let edges_ = mfn.query.Get('_e', vert_);

      if (ifn.len(edges_) == 1) {

        return 0;
      }

      let vecs_ = mfn.calc.Vector(edges_);

      let vec0_ = ifn.vecRev(ifn.vecNorm(vecs_[pythonList(0, vecs_.length)]));

      let vec1_ = ifn.vecNorm(vecs_[pythonList(1, vecs_.length)]);

      return ifn.radToDeg(ifn.vecAng2(vec1_, vec0_, [0, 0, 1]));
    }


    async function exec_cluster1_node_9bpb8dz1lre__touchingEdge_($p, edges_, xyz_) {

      for (let edge_ of edges_) {

        let r_ = mfn.calc.Ray(edge_);

        let d_ = ifn.distance(xyz_, ifn.project(xyz_, r_, 0));

        if (d_ < 1) {

          return edge_;
        }
      }

      return null;
    }


    async function exec_cluster1_node_9bpb8dz1lre__trimLoc_($p, block_, parts_) {

      let plines_loc_ = await exec_cluster1_node_9bpb8dz1lre__getPerpLocPlines_($p, block_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let loc_off_ = mfn.poly2d.OffsetMitre(plines_loc_, PART_LOC_D_, PART_LOC_D_, 'square_end');

      let parts_trim_ = mfn.poly2d.Boolean(parts_, loc_off_, 'difference');

      for (let part_ of parts_) {

        let parts_loc_ = mfn.poly2d.Boolean(loc_off_, part_, 'intersect');

        mfn.list.Add(parts_trim_, parts_loc_, 'extend_end');
      }

      mfn.edit.Delete([plines_loc_, loc_off_, parts_], 'delete_selected');

      return parts_trim_;
    }


    async function exec_cluster1_node_9bpb8dz1lre__getPerimPlines_($p, block_, road_descr_) {

      let posis_ = [];

      for (let edge_ of mfn.query.Get('_e', block_)) {

        let edge_road_ = mfn.attrib.Get(edge_, 'road');

        if (edge_road_ == road_descr_) {

          let start_end_ = mfn.query.Get('ps', edge_);

          if (ifn.len(posis_) == 0 || start_end_[pythonList(0, start_end_.length)] != posis_[pythonList(-1, posis_.length)][pythonList(-1, posis_[pythonList(-1, posis_.length)].length)]) {

            mfn.list.Add(posis_, start_end_, 'to_end');
          } else {

            mfn.list.Add(posis_[pythonList(-1, posis_.length)], start_end_[pythonList(1, start_end_.length)], 'to_end');
          }
        }
      }

      if (ifn.len(posis_) == 0) {

        return [];
      }

      if (ifn.len(posis_) > 1 && posis_[pythonList(-1, posis_.length)][pythonList(-1, posis_[pythonList(-1, posis_.length)].length)] == posis_[pythonList(0, posis_.length)][pythonList(0, posis_[pythonList(0, posis_.length)].length)]) {

        let first_list_ = ifn.listJoin(posis_[pythonList(-1, posis_.length)], posis_[pythonList(0, posis_.length)].slice(1));

        posis_[pythonList(0, posis_.length)] = first_list_;

        posis_ = posis_.slice(0, -1);
      }

      let site_plines_ = mfn.make.Polyline(posis_, 'open');

      await exec_cluster1_node_9bpb8dz1lre__extendPline_($p, site_plines_, 100);
      if ($p.terminated) {
        return mfn.getModel();
      }

      return site_plines_;
    }


    async function exec_cluster1_node_9bpb8dz1lre__getPerpLocPlines_($p, block_) {

      let perp_dir_ = mfn.attrib.Get(null, ['ortho', 2]);

      if (mfn.attrib.Get(block_, 'block_type') == "sec") {

        perp_dir_ = mfn.attrib.Get(null, ['ortho', 1]);
      }

      let loc_perp_edges_ = [];

      for (let edge_ of mfn.query.Get('_e', block_)) {

        let edge_road_ = mfn.attrib.Get(edge_, 'road');

        let edge_vec_ = mfn.calc.Vector(edge_);

        if (edge_road_ == "road_loc" && ifn.abs(ifn.vecDot(ifn.vecNorm(edge_vec_), perp_dir_)) > 0.8) {

          mfn.list.Add(loc_perp_edges_, edge_, 'to_end');
        }
      }

      let plines_ = mfn.make.Polyline(loc_perp_edges_, 'open');

      await exec_cluster1_node_9bpb8dz1lre__extendPline_($p, plines_, 100);
      if ($p.terminated) {
        return mfn.getModel();
      }

      return plines_;
    }


    async function exec_cluster1_node_9bpb8dz1lre__extendPline_($p, plines_, dist_) {

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

    async function exec_cluster1_node_9bpb8dz1lre($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: generate art sec parts no offgrid ', '__null__')
      }


      let all_blocks_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', "block");

      let art_blocks_ = mfn.query.Filter(all_blocks_, ['block_type', null], '==', "art");

      let sec_blocks_ = mfn.query.Filter(all_blocks_, ['block_type', null], '==', "sec");

      mfn.edit.Delete([art_blocks_, sec_blocks_], 'keep_selected');

      let depths_dict_ = {
        "road_art": PART_ART_D_,
        "road_sec": PART_SEC_D_,
        "road_loc": PART_LOC_D_,
        "cold": 0
      };

      let areas_dict_ = {};

      areas_dict_["art_art"] = PART_ART_D_ * PART_ART_D_;

      areas_dict_["art_sec"] = PART_ART_D_ * PART_SEC_D_;

      areas_dict_["art_loc"] = PART_ART_D_ * PART_LOC_D_;

      areas_dict_["sec_sec"] = PART_SEC_D_ * PART_SEC_D_;

      areas_dict_["sec_loc"] = PART_SEC_D_ * PART_LOC_D_;

      for (let block_ of mfn.query.Get('pg', null)) {

        let parts_ = await exec_cluster1_node_9bpb8dz1lre_createParts_($p, block_, depths_dict_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(parts_, `class`, "part");

        await exec_cluster1_node_9bpb8dz1lre_copyAttribs_($p, block_, parts_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        await exec_cluster1_node_9bpb8dz1lre_transferEdgeAttribsTouching_($p, mfn.query.Get('_e', block_), mfn.query.Get('_e', parts_));
        if ($p.terminated) {
          return mfn.getModel();
        }

        for (let part_ of parts_) {

          await exec_cluster1_node_9bpb8dz1lre_addTypeAttribsCrvRoads_($p, part_);
          if ($p.terminated) {
            return mfn.getModel();
          }

          let type_ = mfn.attrib.Get(part_, 'type');

          let corner_expected_area_ = areas_dict_[pythonList(type_, areas_dict_.length)];

          if (corner_expected_area_ != undefined) {

            let area_ = mfn.calc.Area(part_);

            if (area_ > (corner_expected_area_ * 3)) {

              mfn.attrib.Set(part_, `type`, type_.slice(0, 3));
            }
          }
        }

        await exec_cluster1_node_9bpb8dz1lre_transferEdgeAttribsBtwTouchingParts_($p, parts_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.edit.Delete(block_, 'delete_selected');
      }

      for (let part_ of mfn.query.Get('pg', null)) {

        let type_ = mfn.attrib.Get(part_, 'type');

        if (ifn.len(type_) == 7) {

          let lengths_ = mfn.calc.Length(mfn.query.Filter(mfn.query.Get('_e', part_), ['road', null], '==', "road_loc"));

          if (ifn.max(lengths_) > (PART_LOC_D_ * 2)) {

            mfn.attrib.Set(part_, `type`, type_.slice(0, 3));
          }

          lengths_ = mfn.calc.Length(mfn.query.Filter(mfn.query.Get('_e', part_), ['road', null], '==', "road_sec"));

          if (ifn.max(lengths_) > (PART_SEC_D_ * 2)) {

            mfn.attrib.Set(part_, `type`, type_.slice(0, 3));
          }
        }
      }
    }


    async function exec_cluster1_node_1bbah65kkdi_createParts_($p, block_, depths_dict_) {

      let quad_ = await exec_cluster1_node_1bbah65kkdi__convertToQuad_($p, block_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      if (quad_ == null) {

        mfn.visualize.Color(block_, [0, 1, 0]);

        mfn.attrib.Set(block_, `type`, "park");

        return [];
      }

      let blk_edges_ = mfn.query.Get('_e', block_);

      let quad_edges_ = mfn.query.Get('_e', quad_);

      await exec_cluster1_node_1bbah65kkdi__transferEdgeAttribs_($p, blk_edges_, quad_edges_, false);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let roads_ = mfn.attrib.Get(quad_edges_, 'road');

      let lengths_ = mfn.calc.Length(quad_edges_);

      let rays_ = mfn.calc.Ray(quad_edges_);

      let a_ = ifn.listSort(ifn.range(4), lengths_)[pythonList(0, ifn.listSort(ifn.range(4), lengths_).length)];

      let b_ = (a_ + 1) % 4;

      let c_ = (a_ + 2) % 4;

      let d_ = (a_ + 3) % 4;

      let off_a_ = depths_dict_[pythonList(roads_[pythonList(a_, roads_.length)], depths_dict_.length)];

      let off_b_ = depths_dict_[pythonList(roads_[pythonList(b_, roads_.length)], depths_dict_.length)];

      let off_c_ = depths_dict_[pythonList(roads_[pythonList(c_, roads_.length)], depths_dict_.length)];

      let off_d_ = depths_dict_[pythonList(roads_[pythonList(d_, roads_.length)], depths_dict_.length)];

      let thick_ = off_b_ + off_d_;

      let thick_70_ = thick_ * 0.7;

      if (lengths_[pythonList(a_, lengths_.length)] < thick_70_ && lengths_[pythonList(c_, lengths_.length)] < thick_70_) {

        mfn.edit.Delete(quad_, 'delete_selected');

        mfn.attrib.Set(block_, `type`, "loc");

        return [block_];
      } else {
        if (lengths_[pythonList(a_, lengths_.length)] < thick_70_) {
        } else {
          if (lengths_[pythonList(c_, lengths_.length)] < thick_70_) {
          }
        }
      }

      let vperp_a_ = ifn.vecSetLen([-rays_[pythonList(a_, rays_.length)][pythonList(1, rays_[pythonList(a_, rays_.length)].length)][pythonList(1, rays_[pythonList(a_, rays_.length)][pythonList(1, rays_[pythonList(a_, rays_.length)].length)].length)], rays_[pythonList(a_, rays_.length)][pythonList(1, rays_[pythonList(a_, rays_.length)].length)][pythonList(0, rays_[pythonList(a_, rays_.length)][pythonList(1, rays_[pythonList(a_, rays_.length)].length)].length)], 0], off_a_);

      let o_a_ = ifn.vecAdd(rays_[pythonList(a_, rays_.length)][pythonList(0, rays_[pythonList(a_, rays_.length)].length)], vperp_a_);

      let ray_a_inn_ = ifn.rayMake(o_a_, rays_[pythonList(a_, rays_.length)][pythonList(1, rays_[pythonList(a_, rays_.length)].length)]);

      let vperp_c_ = ifn.vecSetLen([-rays_[pythonList(c_, rays_.length)][pythonList(1, rays_[pythonList(c_, rays_.length)].length)][pythonList(1, rays_[pythonList(c_, rays_.length)][pythonList(1, rays_[pythonList(c_, rays_.length)].length)].length)], rays_[pythonList(c_, rays_.length)][pythonList(1, rays_[pythonList(c_, rays_.length)].length)][pythonList(0, rays_[pythonList(c_, rays_.length)][pythonList(1, rays_[pythonList(c_, rays_.length)].length)].length)], 0], off_c_);

      let o_c_ = ifn.vecAdd(rays_[pythonList(c_, rays_.length)][pythonList(0, rays_[pythonList(c_, rays_.length)].length)], vperp_c_);

      let ray_c_inn_ = ifn.rayMake(o_c_, rays_[pythonList(c_, rays_.length)][pythonList(1, rays_[pythonList(c_, rays_.length)].length)]);

      let isect_ab_ = ifn.intersect(rays_[pythonList(b_, rays_.length)], ray_a_inn_, 2);

      let isect_ad_ = ifn.intersect(rays_[pythonList(d_, rays_.length)], ray_a_inn_, 2);

      let isect_cb_ = ifn.intersect(rays_[pythonList(b_, rays_.length)], ray_c_inn_, 2);

      let isect_cd_ = ifn.intersect(rays_[pythonList(d_, rays_.length)], ray_c_inn_, 2);

      let xyz_a_ = ifn.vecDiv(ifn.vecAdd(isect_ab_, isect_ad_), 2);

      let xyz_c_ = ifn.vecDiv(ifn.vecAdd(isect_cb_, isect_cd_), 2);

      let parts_ = [];

      if (ifn.distance(xyz_a_, xyz_c_) < 20) {

        let xyz_cen_ = ifn.vecDiv(ifn.vecAdd(xyz_a_, xyz_c_), 2);

        parts_ = await exec_cluster1_node_1bbah65kkdi_creat4Parts_($p, blk_edges_, quad_edges_, rays_, [a_, b_, c_, d_], xyz_cen_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(parts_, `type`, "loc_loc");
      } else {
        if (lengths_[pythonList(a_, lengths_.length)] < thick_70_ || lengths_[pythonList(c_, lengths_.length)] < thick_70_) {

          parts_ = await exec_cluster1_node_1bbah65kkdi_creat5Parts_($p, blk_edges_, quad_edges_, rays_, [a_, b_, c_, d_], [xyz_a_, xyz_c_], lengths_[pythonList(a_, lengths_.length)] < lengths_[pythonList(c_, lengths_.length)]);
          if ($p.terminated) {
            return mfn.getModel();
          }

          mfn.attrib.Set(parts_.slice(0, 3), `type`, "loc_loc");

          mfn.attrib.Set(parts_.slice(3), `type`, "loc");
        } else {

          parts_ = await exec_cluster1_node_1bbah65kkdi_creat6Parts_($p, blk_edges_, quad_edges_, rays_, [a_, b_, c_, d_], [xyz_a_, xyz_c_]);
          if ($p.terminated) {
            return mfn.getModel();
          }

          mfn.attrib.Set(parts_.slice(0, 4), `type`, "loc_loc");

          mfn.attrib.Set(parts_.slice(4), `type`, "loc");
        }
      }

      if (ifn.len(blk_edges_) != 4) {

        let trimmed_ = mfn.poly2d.Boolean(parts_, block_, 'intersect');

        mfn.edit.Delete(parts_, 'delete_selected');

        parts_ = trimmed_;
      }

      for (let part_ of parts_) {

        await exec_cluster1_node_1bbah65kkdi__transferEdgeAttribs_($p, blk_edges_, mfn.query.Get('_e', part_), true);
        if ($p.terminated) {
          return mfn.getModel();
        }
      }

      await exec_cluster1_node_1bbah65kkdi__copyAttribs_($p, block_, parts_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      mfn.edit.Delete([block_, quad_], 'delete_selected');

      return parts_;
    }


    async function exec_cluster1_node_1bbah65kkdi_creat4Parts_($p, blk_edges_, quad_edges_, rays_, abcd_, xyz_cen_) {

      let a_ = abcd_[pythonList(0, abcd_.length)];

      let b_ = abcd_[pythonList(1, abcd_.length)];

      let c_ = abcd_[pythonList(2, abcd_.length)];

      let d_ = abcd_[pythonList(3, abcd_.length)];

      let xyz_ab_ = ifn.project(xyz_cen_, rays_[pythonList(b_, rays_.length)]);

      let xyz_ad_ = ifn.project(xyz_cen_, rays_[pythonList(d_, rays_.length)]);

      if (ifn.len(blk_edges_) != 4) {

        for (let i_ of ifn.range(4)) {

          let r_o_ = rays_[pythonList(i_, rays_.length)][pythonList(0, rays_[pythonList(i_, rays_.length)].length)];

          let r_v_ = rays_[pythonList(i_, rays_.length)][pythonList(1, rays_[pythonList(i_, rays_.length)].length)];

          let r_vp_ = ifn.vecSetLen([r_v_[pythonList(1, r_v_.length)], -r_v_[pythonList(0, r_v_.length)], 0], 10);

          rays_[pythonList(i_, rays_.length)] = ifn.rayMake(ifn.vecAdd(r_o_, r_vp_), r_v_);
        }
      }

      let a_mid_ = ifn.vecAdd(rays_[pythonList(a_, rays_.length)][pythonList(0, rays_[pythonList(a_, rays_.length)].length)], ifn.vecDiv(rays_[pythonList(a_, rays_.length)][pythonList(1, rays_[pythonList(a_, rays_.length)].length)], 2));

      let c_mid_ = ifn.vecAdd(rays_[pythonList(c_, rays_.length)][pythonList(0, rays_[pythonList(c_, rays_.length)].length)], ifn.vecDiv(rays_[pythonList(c_, rays_.length)][pythonList(1, rays_[pythonList(c_, rays_.length)].length)], 2));

      let xyz_ea0_ = ifn.intersect(rays_[pythonList(d_, rays_.length)], rays_[pythonList(a_, rays_.length)], 2);

      let xyz_ea1_ = ifn.intersect(rays_[pythonList(a_, rays_.length)], rays_[pythonList(b_, rays_.length)], 2);

      let xyz_ec0_ = ifn.intersect(rays_[pythonList(b_, rays_.length)], rays_[pythonList(c_, rays_.length)], 2);

      let xyz_ec1_ = ifn.intersect(rays_[pythonList(c_, rays_.length)], rays_[pythonList(d_, rays_.length)], 2);

      let posis_a_ = mfn.make.Position([xyz_ea0_, xyz_ea1_]);

      let posis_c_ = mfn.make.Position([xyz_ec0_, xyz_ec1_]);

      let posis_cen_ = mfn.make.Position([xyz_ad_, xyz_cen_, xyz_ab_]);

      let posis_a_mid_ = mfn.make.Position(a_mid_);

      let posis_c_mid_ = mfn.make.Position(c_mid_);

      let part_ab_ = mfn.make.Polygon([posis_a_[pythonList(1, posis_a_.length)], posis_cen_[pythonList(2, posis_cen_.length)], posis_cen_[pythonList(1, posis_cen_.length)], posis_a_mid_]);

      let part_ad_ = mfn.make.Polygon([posis_a_[pythonList(0, posis_a_.length)], posis_a_mid_, posis_cen_[pythonList(1, posis_cen_.length)], posis_cen_[pythonList(0, posis_cen_.length)]]);

      let part_cb_ = mfn.make.Polygon([posis_c_[pythonList(1, posis_c_.length)], posis_cen_[pythonList(0, posis_cen_.length)], posis_cen_[pythonList(1, posis_cen_.length)], posis_c_mid_]);

      let part_cd_ = mfn.make.Polygon([posis_c_[pythonList(0, posis_c_.length)], posis_c_mid_, posis_cen_[pythonList(1, posis_cen_.length)], posis_cen_[pythonList(2, posis_cen_.length)]]);

      let parts_ = [part_ab_, part_ad_, part_cb_, part_cd_];

      return parts_;
    }


    async function exec_cluster1_node_1bbah65kkdi_creat5Parts_($p, blk_edges_, quad_edges_, rays_, abcd_, xyzs_, a_is_short_) {

      let a_ = abcd_[pythonList(0, abcd_.length)];

      let b_ = abcd_[pythonList(1, abcd_.length)];

      let c_ = abcd_[pythonList(2, abcd_.length)];

      let d_ = abcd_[pythonList(3, abcd_.length)];

      let xyz_a_ = xyzs_[pythonList(0, xyzs_.length)];

      let xyz_c_ = xyzs_[pythonList(1, xyzs_.length)];

      if (ifn.len(blk_edges_) != 4) {

        for (let i_ of ifn.range(4)) {

          let r_o_ = rays_[pythonList(i_, rays_.length)][pythonList(0, rays_[pythonList(i_, rays_.length)].length)];

          let r_v_ = rays_[pythonList(i_, rays_.length)][pythonList(1, rays_[pythonList(i_, rays_.length)].length)];

          let r_vp_ = ifn.vecSetLen([r_v_[pythonList(1, r_v_.length)], -r_v_[pythonList(0, r_v_.length)], 0], 10);

          rays_[pythonList(i_, rays_.length)] = ifn.rayMake(ifn.vecAdd(r_o_, r_vp_), r_v_);
        }
      }

      let xyz_ab_ = ifn.project(xyz_a_, rays_[pythonList(b_, rays_.length)]);

      let xyz_ad_ = ifn.project(xyz_a_, rays_[pythonList(d_, rays_.length)]);

      let xyz_cb_ = ifn.project(xyz_c_, rays_[pythonList(b_, rays_.length)]);

      let xyz_cd_ = ifn.project(xyz_c_, rays_[pythonList(d_, rays_.length)]);

      let posis_a_ = mfn.query.Get('ps', quad_edges_[pythonList(a_, quad_edges_.length)]);

      let posis_c_ = mfn.query.Get('ps', quad_edges_[pythonList(c_, quad_edges_.length)]);

      let xyz_ea0_ = ifn.intersect(rays_[pythonList(d_, rays_.length)], rays_[pythonList(a_, rays_.length)], 2);

      let xyz_ea1_ = ifn.intersect(rays_[pythonList(a_, rays_.length)], rays_[pythonList(b_, rays_.length)], 2);

      let xyz_ec0_ = ifn.intersect(rays_[pythonList(b_, rays_.length)], rays_[pythonList(c_, rays_.length)], 2);

      let xyz_ec1_ = ifn.intersect(rays_[pythonList(c_, rays_.length)], rays_[pythonList(d_, rays_.length)], 2);

      posis_a_ = mfn.make.Position([xyz_ea0_, xyz_ea1_]);

      posis_c_ = mfn.make.Position([xyz_ec0_, xyz_ec1_]);

      let posis_off_a_ = mfn.make.Position([xyz_ad_, xyz_a_, xyz_ab_]);

      let posis_off_c_ = mfn.make.Position([xyz_cb_, xyz_c_, xyz_cd_]);

      let parts_ = [];

      if (a_is_short_) {

        let c_mid_ = ifn.vecAdd(rays_[pythonList(c_, rays_.length)][pythonList(0, rays_[pythonList(c_, rays_.length)].length)], ifn.vecDiv(rays_[pythonList(c_, rays_.length)][pythonList(1, rays_[pythonList(c_, rays_.length)].length)], 2));

        let posis_c_mid_ = mfn.make.Position(c_mid_);

        let part_a_ = mfn.make.Polygon([posis_a_[pythonList(1, posis_a_.length)], posis_off_a_[pythonList(2, posis_off_a_.length)], posis_off_a_[pythonList(1, posis_off_a_.length)], posis_off_a_[pythonList(0, posis_off_a_.length)], posis_a_[pythonList(0, posis_a_.length)]]);

        let part_cb_ = mfn.make.Polygon([posis_c_[pythonList(1, posis_c_.length)], posis_off_c_[pythonList(2, posis_off_c_.length)], posis_off_c_[pythonList(1, posis_off_c_.length)], posis_c_mid_]);

        let part_cd_ = mfn.make.Polygon([posis_c_[pythonList(0, posis_c_.length)], posis_c_mid_, posis_off_c_[pythonList(1, posis_off_c_.length)], posis_off_c_[pythonList(0, posis_off_c_.length)]]);

        let part_b_ = mfn.make.Polygon([posis_off_a_[pythonList(2, posis_off_a_.length)], posis_off_c_[pythonList(0, posis_off_c_.length)], posis_off_c_[pythonList(1, posis_off_c_.length)], posis_off_a_[pythonList(1, posis_off_a_.length)]]);

        let part_d_ = mfn.make.Polygon([posis_off_c_[pythonList(2, posis_off_c_.length)], posis_off_a_[pythonList(0, posis_off_a_.length)], posis_off_a_[pythonList(1, posis_off_a_.length)], posis_off_c_[pythonList(1, posis_off_c_.length)]]);

        parts_ = [part_a_, part_cb_, part_cd_, part_b_, part_d_];
      } else {

        let a_mid_ = ifn.vecAdd(rays_[pythonList(a_, rays_.length)][pythonList(0, rays_[pythonList(a_, rays_.length)].length)], ifn.vecDiv(rays_[pythonList(a_, rays_.length)][pythonList(1, rays_[pythonList(a_, rays_.length)].length)], 2));

        let posis_a_mid_ = mfn.make.Position(a_mid_);

        let part_ab_ = mfn.make.Polygon([posis_a_[pythonList(1, posis_a_.length)], posis_off_a_[pythonList(2, posis_off_a_.length)], posis_off_a_[pythonList(1, posis_off_a_.length)], posis_a_mid_]);

        let part_ad_ = mfn.make.Polygon([posis_a_[pythonList(0, posis_a_.length)], posis_a_mid_, posis_off_a_[pythonList(1, posis_off_a_.length)], posis_off_a_[pythonList(0, posis_off_a_.length)]]);

        let part_c_ = mfn.make.Polygon([posis_c_[pythonList(1, posis_c_.length)], posis_off_c_[pythonList(2, posis_off_c_.length)], posis_off_c_[pythonList(1, posis_off_c_.length)], posis_off_c_[pythonList(0, posis_off_c_.length)], posis_c_[pythonList(0, posis_c_.length)]]);

        let part_b_ = mfn.make.Polygon([posis_off_a_[pythonList(2, posis_off_a_.length)], posis_off_c_[pythonList(0, posis_off_c_.length)], posis_off_c_[pythonList(1, posis_off_c_.length)], posis_off_a_[pythonList(1, posis_off_a_.length)]]);

        let part_d_ = mfn.make.Polygon([posis_off_c_[pythonList(2, posis_off_c_.length)], posis_off_a_[pythonList(0, posis_off_a_.length)], posis_off_a_[pythonList(1, posis_off_a_.length)], posis_off_c_[pythonList(1, posis_off_c_.length)]]);

        parts_ = [part_ab_, part_ad_, part_c_, part_b_, part_d_];
      }

      return parts_;
    }


    async function exec_cluster1_node_1bbah65kkdi_creat6Parts_($p, blk_edges_, quad_edges_, rays_, abcd_, xyzs_) {

      let a_ = abcd_[pythonList(0, abcd_.length)];

      let b_ = abcd_[pythonList(1, abcd_.length)];

      let c_ = abcd_[pythonList(2, abcd_.length)];

      let d_ = abcd_[pythonList(3, abcd_.length)];

      let xyz_a_ = xyzs_[pythonList(0, xyzs_.length)];

      let xyz_c_ = xyzs_[pythonList(1, xyzs_.length)];

      if (ifn.len(blk_edges_) != 4) {

        for (let i_ of ifn.range(4)) {

          let r_o_ = rays_[pythonList(i_, rays_.length)][pythonList(0, rays_[pythonList(i_, rays_.length)].length)];

          let r_v_ = rays_[pythonList(i_, rays_.length)][pythonList(1, rays_[pythonList(i_, rays_.length)].length)];

          let r_vp_ = ifn.vecSetLen([r_v_[pythonList(1, r_v_.length)], -r_v_[pythonList(0, r_v_.length)], 0], 10);

          rays_[pythonList(i_, rays_.length)] = ifn.rayMake(ifn.vecAdd(r_o_, r_vp_), r_v_);
        }
      }

      let xyz_ab_ = ifn.project(xyz_a_, rays_[pythonList(b_, rays_.length)]);

      let xyz_ad_ = ifn.project(xyz_a_, rays_[pythonList(d_, rays_.length)]);

      let xyz_cb_ = ifn.project(xyz_c_, rays_[pythonList(b_, rays_.length)]);

      let xyz_cd_ = ifn.project(xyz_c_, rays_[pythonList(d_, rays_.length)]);

      let a_mid_ = ifn.vecAdd(rays_[pythonList(a_, rays_.length)][pythonList(0, rays_[pythonList(a_, rays_.length)].length)], ifn.vecDiv(rays_[pythonList(a_, rays_.length)][pythonList(1, rays_[pythonList(a_, rays_.length)].length)], 2));

      let c_mid_ = ifn.vecAdd(rays_[pythonList(c_, rays_.length)][pythonList(0, rays_[pythonList(c_, rays_.length)].length)], ifn.vecDiv(rays_[pythonList(c_, rays_.length)][pythonList(1, rays_[pythonList(c_, rays_.length)].length)], 2));

      let xyz_ea0_ = ifn.intersect(rays_[pythonList(d_, rays_.length)], rays_[pythonList(a_, rays_.length)], 2);

      let xyz_ea1_ = ifn.intersect(rays_[pythonList(a_, rays_.length)], rays_[pythonList(b_, rays_.length)], 2);

      let xyz_ec0_ = ifn.intersect(rays_[pythonList(b_, rays_.length)], rays_[pythonList(c_, rays_.length)], 2);

      let xyz_ec1_ = ifn.intersect(rays_[pythonList(c_, rays_.length)], rays_[pythonList(d_, rays_.length)], 2);

      let posis_a_ = mfn.make.Position([xyz_ea0_, xyz_ea1_]);

      let posis_c_ = mfn.make.Position([xyz_ec0_, xyz_ec1_]);

      let posis_off_a_ = mfn.make.Position([xyz_ad_, xyz_a_, xyz_ab_]);

      let posis_off_c_ = mfn.make.Position([xyz_cb_, xyz_c_, xyz_cd_]);

      let posis_a_mid_ = mfn.make.Position(a_mid_);

      let posis_c_mid_ = mfn.make.Position(c_mid_);

      let part_ab_ = mfn.make.Polygon([posis_a_[pythonList(1, posis_a_.length)], posis_off_a_[pythonList(2, posis_off_a_.length)], posis_off_a_[pythonList(1, posis_off_a_.length)], posis_a_mid_]);

      let part_ad_ = mfn.make.Polygon([posis_a_[pythonList(0, posis_a_.length)], posis_a_mid_, posis_off_a_[pythonList(1, posis_off_a_.length)], posis_off_a_[pythonList(0, posis_off_a_.length)]]);

      let part_cb_ = mfn.make.Polygon([posis_c_[pythonList(1, posis_c_.length)], posis_off_c_[pythonList(2, posis_off_c_.length)], posis_off_c_[pythonList(1, posis_off_c_.length)], posis_c_mid_]);

      let part_cd_ = mfn.make.Polygon([posis_c_[pythonList(0, posis_c_.length)], posis_c_mid_, posis_off_c_[pythonList(1, posis_off_c_.length)], posis_off_c_[pythonList(0, posis_off_c_.length)]]);

      let part_b_ = mfn.make.Polygon([posis_off_a_[pythonList(2, posis_off_a_.length)], posis_off_c_[pythonList(0, posis_off_c_.length)], posis_off_c_[pythonList(1, posis_off_c_.length)], posis_off_a_[pythonList(1, posis_off_a_.length)]]);

      let part_d_ = mfn.make.Polygon([posis_off_c_[pythonList(2, posis_off_c_.length)], posis_off_a_[pythonList(0, posis_off_a_.length)], posis_off_a_[pythonList(1, posis_off_a_.length)], posis_off_c_[pythonList(1, posis_off_c_.length)]]);

      let parts_ = [part_ab_, part_ad_, part_cb_, part_cd_, part_b_, part_d_];

      return parts_;
    }


    async function exec_cluster1_node_1bbah65kkdi_addTypeAttribsCrvRoads_($p, parts_) {

      for (let part_ of parts_) {

        let part_edges_ = mfn.query.Get('_e', part_);

        let road_types_ = [];

        let art_ = await exec_cluster1_node_1bbah65kkdi__getCrvRoadEdgeTypes_($p, part_edges_, "road_art");
        if ($p.terminated) {
          return mfn.getModel();
        }

        let sec_ = await exec_cluster1_node_1bbah65kkdi__getCrvRoadEdgeTypes_($p, part_edges_, "road_sec");
        if ($p.terminated) {
          return mfn.getModel();
        }

        let loc_ = await exec_cluster1_node_1bbah65kkdi__getCrvRoadEdgeTypes_($p, part_edges_, "road_loc");
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.list.Add(road_types_, ifn.listFlat([art_, sec_, loc_]), 'extend_end');

        let check_ = road_types_;

        mfn.attrib.Set(part_, `class`, "part");

        let block_types_dict_ = {
          "road_art": "art",
          "road_sec": "sec",
          "road_ter": "ter",
          "road_loc": "loc"
        };

        let cats_dict_ = {
          "road_art": 1,
          "road_sec": 2,
          "road_ter": 3,
          "road_loc": 4
        };

        if (ifn.len(road_types_) == 0) {

          mfn.attrib.Set(part_, `type`, "off_grid");
        } else {
          if (ifn.len(road_types_) == 1) {

            let road_type_ = road_types_[pythonList(0, road_types_.length)];

            let block_type_ = ifn.string(block_types_dict_[pythonList(road_type_, block_types_dict_.length)]);

            mfn.attrib.Set(part_, `type`, block_type_);
          } else {

            let road_type0_ = road_types_[pythonList(0, road_types_.length)];

            let block_type0_ = ifn.string(block_types_dict_[pythonList(road_type0_, block_types_dict_.length)]);

            let road_type1_ = road_types_[pythonList(1, road_types_.length)];

            let block_type1_ = ifn.string(block_types_dict_[pythonList(road_type1_, block_types_dict_.length)]);

            mfn.attrib.Set(part_, `type`, block_type0_ + "_" + block_type1_);
          }
        }
      }
    }


    async function exec_cluster1_node_1bbah65kkdi_transferEdgeAttribsBtwTouchingParts_($p, parts_) {

      let edges_ = mfn.query.Get('_e', parts_);

      for (let to_edge_ of edges_) {

        if (mfn.attrib.Get(to_edge_, 'road') == undefined) {

          let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

          let idx_ = ifn.listFind(edges_, to_edge_);

          let from_edges_ = ifn.listJoin(edges_.slice(0, idx_), edges_.slice(idx_ + 1));

          let from_edge_ = await exec_cluster1_node_1bbah65kkdi__touchingEdge_($p, from_edges_, cen_);
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


    async function exec_cluster1_node_1bbah65kkdi__copyAttribs_($p, pgon_from_, pgons_to_) {

      let block_id_ = mfn.attrib.Get(pgon_from_, 'block_id');

      let block_type_ = mfn.attrib.Get(pgon_from_, 'block_type');

      let site_ = mfn.attrib.Get(pgon_from_, 'site');

      for (let pgon_to_ of pgons_to_) {

        mfn.attrib.Set(pgon_to_, `block_id`, block_id_);

        mfn.attrib.Set(pgon_to_, `block_type`, block_type_);

        mfn.attrib.Set(pgon_to_, `site`, site_);
      }
    }


    async function exec_cluster1_node_1bbah65kkdi__getCrvRoadEdgeTypes_($p, part_edges_, road_descr_) {

      let edges_ = [[]];

      for (let edge_ of part_edges_) {

        let edge_road_ = mfn.attrib.Get(edge_, 'road');

        if (edge_road_ == road_descr_) {

          mfn.list.Add(edges_[pythonList(-1, edges_.length)], edge_, 'to_end');
        } else {

          if (ifn.len(edges_[pythonList(-1, edges_.length)]) != 0) {

            mfn.list.Add(edges_, [], 'to_end');
          }
        }
      }

      if (ifn.len(edges_) == 2) {

        edges_ = ifn.listJoin(edges_[pythonList(1, edges_.length)], edges_[pythonList(0, edges_.length)]);
      } else {

        edges_ = edges_[pythonList(0, edges_.length)];
      }

      if (ifn.len(edges_) == 0) {

        return [];
      }

      if (ifn.len(edges_) == 1) {

        return [road_descr_];
      }

      let road_types_ = [road_descr_];

      for (let i_ of ifn.range(1, ifn.len(edges_))) {

        let ang_ = await exec_cluster1_node_1bbah65kkdi__vertAng_($p, mfn.query.Get('_v', edges_[pythonList(i_, edges_.length)])[pythonList(0, mfn.query.Get('_v', edges_[pythonList(i_, edges_.length)]).length)]);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (ifn.abs(ang_ - 180) > 45) {

          mfn.list.Add(road_types_, road_descr_, 'to_end');
        }
      }

      return road_types_;
    }


    async function exec_cluster1_node_1bbah65kkdi__convertToQuad_($p, pgon_) {

      let area_ = mfn.calc.Area(pgon_);

      if (area_ < 100) {

        return null;
      } else {
        if (ifn.len(mfn.query.Get('_v', pgon_)) < 4) {

          return null;
        }
      }

      let pgon2_ = mfn.make.Copy(pgon_, null);

      let verts_ = mfn.query.Get('_v', pgon2_);

      if (ifn.len(verts_) == 4) {

        return pgon2_;
      }

      let dots_ = [];

      for (let vert_ of verts_) {

        let dot_ = await exec_cluster1_node_1bbah65kkdi__angDot_($p, vert_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.list.Add(dots_, dot_, 'to_end');
      }

      let sorted_verts_ = ifn.listRev(ifn.listSort(ifn.listZip([verts_, dots_]), dots_));

      if (ifn.abs(sorted_verts_[pythonList(ifn.len(verts_) - 4, sorted_verts_.length)][pythonList(1, sorted_verts_[pythonList(ifn.len(verts_) - 4, sorted_verts_.length)].length)]) > 0.99) {

        mfn.edit.Delete(pgon2_, 'delete_selected');

        return null;
      }

      for (let i_ of ifn.range(ifn.len(verts_) - 4)) {

        mfn.edit.Delete(mfn.query.Get('ps', sorted_verts_[pythonList(i_, sorted_verts_.length)][pythonList(0, sorted_verts_[pythonList(i_, sorted_verts_.length)].length)]), 'delete_selected');
      }

      return pgon2_;
    }


    async function exec_cluster1_node_1bbah65kkdi__angDot_($p, vert_) {

      let edges_ = mfn.query.Get('_e', vert_);

      let vecs_ = mfn.calc.Vector(edges_);

      let vec0_ = ifn.vecNorm(vecs_[pythonList(0, vecs_.length)]);

      let vec1_ = ifn.vecNorm(vecs_[pythonList(1, vecs_.length)]);

      return ifn.vecDot(vec0_, vec1_);
    }


    async function exec_cluster1_node_1bbah65kkdi__vertAng_($p, vert_) {

      let edges_ = mfn.query.Get('_e', vert_);

      if (ifn.len(edges_) == 1) {

        return 0;
      }

      let vecs_ = mfn.calc.Vector(edges_);

      let vec0_ = ifn.vecRev(ifn.vecNorm(vecs_[pythonList(0, vecs_.length)]));

      let vec1_ = ifn.vecNorm(vecs_[pythonList(1, vecs_.length)]);

      return ifn.radToDeg(ifn.vecAng2(vec1_, vec0_, [0, 0, 1]));
    }


    async function exec_cluster1_node_1bbah65kkdi__transferEdgeAttribs_($p, from_edges_, to_edges_, must_touch_) {

      for (let to_edge_ of to_edges_) {

        let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

        let from_edge_ = null;

        if (must_touch_) {

          from_edge_ = await exec_cluster1_node_1bbah65kkdi__touchingEdge_($p, from_edges_, cen_);
          if ($p.terminated) {
            return mfn.getModel();
          }
        } else {

          from_edge_ = await exec_cluster1_node_1bbah65kkdi__closestEdge_($p, from_edges_, cen_);
          if ($p.terminated) {
            return mfn.getModel();
          }
        }

        if (from_edge_ != null) {

          mfn.attrib.Set(to_edge_, `road`, mfn.attrib.Get(from_edge_, 'road'));
        }
      }
    }


    async function exec_cluster1_node_1bbah65kkdi__touchingEdge_($p, edges_, xyz_) {

      for (let edge_ of edges_) {

        let r_ = mfn.calc.Ray(edge_);

        let d_ = ifn.distance(xyz_, ifn.project(xyz_, r_, 0));

        if (d_ < 1) {

          return edge_;
        }
      }

      return null;
    }


    async function exec_cluster1_node_1bbah65kkdi__closestEdge_($p, edges_, xyz_) {

      let min_dist_ = Infinity;

      let close_edge_ = null;

      for (let edge_ of edges_) {

        let r_ = mfn.calc.Ray(edge_);

        let d_ = ifn.distance(xyz_, ifn.project(xyz_, r_, 0));

        if (d_ < min_dist_) {

          min_dist_ = d_;

          close_edge_ = edge_;
        }
      }

      return close_edge_;
    }


    async function exec_cluster1_node_1bbah65kkdi__projectOntoEdge_($p, xyz_, edges_) {

      let d_min_ = Infinity;

      let xyz_proj_min_ = xyz_;

      for (let edge_ of edges_) {

        let r_ = mfn.calc.Ray(edge_);

        let xyz_proj_ = ifn.project(xyz_, r_, 0);

        let d_ = ifn.distance(xyz_, xyz_proj_);

        if (d_ < d_min_) {

          d_min_ = d_;

          xyz_proj_min_ = xyz_proj_;
        }
      }

      return xyz_proj_min_;
    }

    async function exec_cluster1_node_1bbah65kkdi($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: generate loc parts no offgrid', '__null__')
      }


      let all_blocks_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', "block");

      let loc_blocks_ = mfn.query.Filter(mfn.query.Get('pg', all_blocks_), ['block_type', null], '==', "loc");

      mfn.edit.Delete(loc_blocks_, 'keep_selected');

      let depths_dict_ = {
        "road_art": PART_ART_D_,
        "road_sec": PART_SEC_D_,
        "road_loc": PART_LOC_D_,
        "cold": 0
      };

      for (let block_ of mfn.query.Get('pg', null)) {

        let parts_ = await exec_cluster1_node_1bbah65kkdi_createParts_($p, block_, depths_dict_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(parts_, `class`, "part");

        await exec_cluster1_node_1bbah65kkdi_addTypeAttribsCrvRoads_($p, parts_);
        if ($p.terminated) {
          return mfn.getModel();
        }
      }
    }


    async function exec_cluster1_node_gi0wr2wsff_subdivOfGrid_($p, off_grid_) {

      let block_id_ = mfn.attrib.Get(off_grid_, 'block_id');

      let block_type_ = mfn.attrib.Get(off_grid_, 'block_type');

      let site_ = mfn.attrib.Get(off_grid_, 'site');

      let swap_orientation_ = mfn.attrib.Get(off_grid_, 'block_type') == "sec";

      let quad_tri_ = await exec_cluster1_node_gi0wr2wsff__convertOffGridToQuadTri_($p, off_grid_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      if (quad_tri_ == null) {

        let copy_ = mfn.make.Copy(off_grid_, null);

        mfn.attrib.Set(copy_, `class`, "plot");

        return [copy_];
      }

      let off_grid_parts_ = off_grid_;

      if (ifn.len(quad_tri_) > 1) {

        off_grid_parts_ = mfn.poly2d.Boolean(off_grid_, quad_tri_[pythonList(1, quad_tri_.length)], 'difference');

        let tri_ = mfn.poly2d.Boolean(quad_tri_[pythonList(1, quad_tri_.length)], off_grid_, 'intersect');

        await exec_cluster1_node_gi0wr2wsff__transferEdgeAttribsTouching_($p, mfn.query.Get('_e', off_grid_), mfn.query.Get('_e', tri_));
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.edit.Delete(quad_tri_[pythonList(1, quad_tri_.length)], 'delete_selected');

        mfn.attrib.Set(tri_, `site`, site_);

        mfn.attrib.Set(tri_, `block_type`, block_type_);

        mfn.attrib.Set(tri_, `block_id`, block_id_);

        mfn.attrib.Set(tri_, `class`, "part");

        mfn.attrib.Set(tri_, `type`, "white");
      }

      let posis_ = await exec_cluster1_node_gi0wr2wsff__createClusterPosis_($p, quad_tri_[pythonList(0, quad_tri_.length)], swap_orientation_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let parts_ = await exec_cluster1_node_gi0wr2wsff__createClusterParts_($p, off_grid_, off_grid_parts_, posis_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      for (let part_ of parts_) {

        mfn.attrib.Set(part_, `site`, site_);

        mfn.attrib.Set(part_, `block_type`, block_type_);

        mfn.attrib.Set(part_, `block_id`, block_id_);

        let area_ = mfn.calc.Area(part_);

        if (area_ > (PART_OG_D_ * PART_OG_W_ * 0.3)) {

          mfn.attrib.Set(part_, `class`, "part");

          mfn.attrib.Set(part_, `type`, "off_grid0");
        } else {

          mfn.attrib.Set(part_, `class`, "park");

          mfn.attrib.Set(part_, `type`, "park");
        }
      }

      mfn.edit.Delete([off_grid_, off_grid_parts_], 'delete_selected');

      await exec_cluster1_node_gi0wr2wsff__transferEdgeAttribsBtwTouchingParts_($p, parts_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      return parts_;
    }


    async function exec_cluster1_node_gi0wr2wsff__convertOffGridToQuadTri_($p, off_grid_) {

      let off_grid_edges_ = mfn.query.Get('_e', off_grid_);

      let area_ = mfn.calc.Area(off_grid_);

      if (area_ < 100) {

        return null;
      } else {
        if (ifn.len(mfn.query.Get('_v', off_grid_)) < 4) {

          return null;
        }
      }

      let pgon2_ = mfn.make.Copy(off_grid_, null);

      let verts_ = mfn.query.Get('_v', pgon2_);

      let dots_ = [];

      for (let vert_ of verts_) {

        let dot_ = await exec_cluster1_node_gi0wr2wsff__angDot_($p, vert_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.list.Add(dots_, dot_, 'to_end');
      }

      let sorted_verts_ = ifn.listRev(ifn.listSort(ifn.listZip([verts_, dots_]), dots_));

      for (let i_ of ifn.range(ifn.len(verts_) - 4)) {

        mfn.edit.Delete(mfn.query.Get('ps', sorted_verts_[pythonList(i_, sorted_verts_.length)][pythonList(0, sorted_verts_[pythonList(i_, sorted_verts_.length)].length)]), 'delete_selected');
      }

      await exec_cluster1_node_gi0wr2wsff__transferEdgeAttribsClosest_($p, off_grid_edges_, mfn.query.Get('_e', pgon2_));
      if ($p.terminated) {
        return mfn.getModel();
      }

      if (sorted_verts_[pythonList(-1, sorted_verts_.length)][pythonList(1, sorted_verts_[pythonList(-1, sorted_verts_.length)].length)] > -0.7) {

        return [pgon2_];
      }

      let worst_vert_ = sorted_verts_[pythonList(-1, sorted_verts_.length)][pythonList(0, sorted_verts_[pythonList(-1, sorted_verts_.length)].length)];

      let worst_posi_ = mfn.query.Get('ps', worst_vert_)[pythonList(0, mfn.query.Get('ps', worst_vert_).length)];

      let edges_ = mfn.query.Get('_e', worst_vert_);

      let rays_ = mfn.calc.Ray(edges_);

      let vec0_ = ifn.vecRev(rays_[pythonList(0, rays_.length)][pythonList(1, rays_[pythonList(0, rays_.length)].length)]);

      let vec1_ = rays_[pythonList(1, rays_.length)][pythonList(1, rays_[pythonList(1, rays_.length)].length)];

      let xyz_ = rays_[pythonList(1, rays_.length)][pythonList(0, rays_[pythonList(1, rays_.length)].length)];

      let new_posi_ = null;

      let chopped_ = null;

      if (ifn.vecLen(vec0_) < ifn.vecLen(vec1_)) {

        let dist_ = ifn.vecDot(vec0_, ifn.vecNorm(vec1_));

        let new_xyz_ = ifn.vecAdd(xyz_, ifn.vecSetLen(vec1_, dist_));

        new_posi_ = mfn.make.Position(new_xyz_);

        chopped_ = await exec_cluster1_node_gi0wr2wsff_makeChopped_($p, new_xyz_, mfn.attrib.Get(mfn.query.Get('ps', edges_[pythonList(0, edges_.length)])[pythonList(0, mfn.query.Get('ps', edges_[pythonList(0, edges_.length)]).length)], 'xyz'));
        if ($p.terminated) {
          return mfn.getModel();
        }
      } else {

        let dist_ = ifn.vecDot(vec1_, ifn.vecNorm(vec0_));

        let new_xyz_ = ifn.vecAdd(xyz_, ifn.vecSetLen(vec0_, dist_));

        new_posi_ = mfn.make.Position(new_xyz_);

        chopped_ = await exec_cluster1_node_gi0wr2wsff_makeChopped_($p, mfn.attrib.Get(mfn.query.Get('ps', edges_[pythonList(1, edges_.length)])[pythonList(1, mfn.query.Get('ps', edges_[pythonList(1, edges_.length)]).length)], 'xyz'), new_xyz_);
        if ($p.terminated) {
          return mfn.getModel();
        }
      }

      let posis_list_ = mfn.query.Get('ps', pgon2_);

      let idx_ = ifn.listFind(posis_list_, worst_posi_);

      posis_list_[pythonList(idx_, posis_list_.length)] = new_posi_;

      let pgon3_ = mfn.make.Polygon(posis_list_);

      mfn.edit.Delete(pgon2_, 'delete_selected');

      return [pgon3_, chopped_];
    }


    async function exec_cluster1_node_gi0wr2wsff_makeChopped_($p, xyz0_, xyz1_) {

      let vec_ = ifn.vecSetLen(ifn.vecFromTo(xyz0_, xyz1_), 100);

      let vec_perp_ = ifn.vecSetLen([-vec_[pythonList(1, vec_.length)], vec_[pythonList(0, vec_.length)], 0], 200);

      let a_ = ifn.vecSub(xyz0_, vec_);

      let b_ = ifn.vecAdd(xyz1_, vec_);

      let c_ = ifn.vecAdd(b_, vec_perp_);

      let d_ = ifn.vecAdd(a_, vec_perp_);

      let posis_ = mfn.make.Position([a_, b_, c_, d_]);

      let pgon_ = mfn.make.Polygon(posis_);

      return pgon_;
    }


    async function exec_cluster1_node_gi0wr2wsff__convertPartToQuad_($p, part_) {

      let off_grid_edges_ = mfn.query.Get('_e', part_);

      let area_ = mfn.calc.Area(part_);

      if (area_ < 100) {

        return null;
      } else {
        if (ifn.len(mfn.query.Get('_v', part_)) < 4) {

          return null;
        }
      }

      let pgon2_ = mfn.make.Copy(part_, null);

      let verts_ = mfn.query.Get('_v', pgon2_);

      if (ifn.len(verts_) == 4) {

        await exec_cluster1_node_gi0wr2wsff__transferEdgeAttribsClosest_($p, off_grid_edges_, mfn.query.Get('_e', pgon2_));
        if ($p.terminated) {
          return mfn.getModel();
        }

        return pgon2_;
      }

      let dots_ = [];

      for (let vert_ of verts_) {

        let dot_ = await exec_cluster1_node_gi0wr2wsff__angDot_($p, vert_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.list.Add(dots_, dot_, 'to_end');
      }

      let sorted_verts_ = ifn.listRev(ifn.listSort(ifn.listZip([verts_, dots_]), dots_));

      for (let i_ of ifn.range(ifn.len(verts_) - 4)) {

        mfn.edit.Delete(mfn.query.Get('ps', sorted_verts_[pythonList(i_, sorted_verts_.length)][pythonList(0, sorted_verts_[pythonList(i_, sorted_verts_.length)].length)]), 'delete_selected');
      }

      await exec_cluster1_node_gi0wr2wsff__transferEdgeAttribsClosest_($p, off_grid_edges_, mfn.query.Get('_e', pgon2_));
      if ($p.terminated) {
        return mfn.getModel();
      }

      return pgon2_;
    }


    async function exec_cluster1_node_gi0wr2wsff__createClusterPosis_($p, quad_, swap_orientation_) {

      let quad_edges_ = mfn.query.Get('_e', quad_);

      let path_plines_ = [];

      let lengths_ = mfn.calc.Length(quad_edges_);

      let rays_ = mfn.calc.Ray(quad_edges_);

      let rays_loc_ = ifn.rayLtoG(rays_, mfn.attrib.Get(null, 'ortho'));

      let vecs_ = [rays_loc_[pythonList(0, rays_loc_.length)][pythonList(1, rays_loc_[pythonList(0, rays_loc_.length)].length)], rays_loc_[pythonList(1, rays_loc_.length)][pythonList(1, rays_loc_[pythonList(1, rays_loc_.length)].length)], rays_loc_[pythonList(2, rays_loc_.length)][pythonList(1, rays_loc_[pythonList(2, rays_loc_.length)].length)], rays_loc_[pythonList(3, rays_loc_.length)][pythonList(1, rays_loc_[pythonList(3, rays_loc_.length)].length)]];

      let vecs_n_ = ifn.vecNorm(vecs_);

      let dots_ = ifn.vecDot([1, 0, 0], vecs_n_);

      let sorted_ = ifn.listSort(ifn.range(4), dots_);

      let a_ = sorted_[pythonList(-1, sorted_.length)];

      let b_ = (a_ + 1) % 4;

      let c_ = (a_ + 2) % 4;

      let d_ = (a_ + 3) % 4;

      rays_ = mfn.calc.Ray(quad_edges_);

      let p0_ = rays_[pythonList(a_, rays_.length)][pythonList(0, rays_[pythonList(a_, rays_.length)].length)];

      let p1_ = rays_[pythonList(b_, rays_.length)][pythonList(0, rays_[pythonList(b_, rays_.length)].length)];

      let p2_ = rays_[pythonList(c_, rays_.length)][pythonList(0, rays_[pythonList(c_, rays_.length)].length)];

      let p3_ = rays_[pythonList(d_, rays_.length)][pythonList(0, rays_[pythonList(d_, rays_.length)].length)];

      let wd_ = [PART_OG_W_, PART_OG_D_];

      if (swap_orientation_) {

        wd_ = [PART_OG_D_, PART_OG_W_];
      }

      let num_x_ = ifn.floor((3 + ifn.min([lengths_[pythonList(a_, lengths_.length)], lengths_[pythonList(c_, lengths_.length)]])) / wd_[pythonList(0, wd_.length)]) + 1;

      let num_y_ = ifn.floor((3 + ifn.min([lengths_[pythonList(b_, lengths_.length)], lengths_[pythonList(d_, lengths_.length)]])) / wd_[pythonList(1, wd_.length)]) + 1;

      if (num_x_ == 1) {

        num_x_ = 2;
      }

      if (num_y_ == 1) {

        num_y_ = 2;
      }

      let vec_right0_ = ifn.vecDiv(rays_[pythonList(a_, rays_.length)][pythonList(1, rays_[pythonList(a_, rays_.length)].length)], num_x_ - 1);

      let vec_right1_ = ifn.vecRev(ifn.vecDiv(rays_[pythonList(c_, rays_.length)][pythonList(1, rays_[pythonList(c_, rays_.length)].length)], num_x_ - 1));

      let vec_right_exp_ = ifn.vecSetLen(vec_right1_, 100);

      let posis_ = [];

      let cluster_parts_ = [];

      for (let i_ of ifn.range(num_x_)) {

        let posis_list_ = [];

        let x0_ = ifn.vecAdd(p0_, ifn.vecMult(vec_right0_, i_));

        let x1_ = ifn.vecAdd(p3_, ifn.vecMult(vec_right1_, i_));

        if (i_ == 0) {

          x0_ = ifn.vecSub(x0_, vec_right_exp_);

          x1_ = ifn.vecSub(x1_, vec_right_exp_);
        } else {
          if (i_ == num_x_ - 1) {

            x0_ = ifn.vecAdd(x0_, vec_right_exp_);

            x1_ = ifn.vecAdd(x1_, vec_right_exp_);
          }
        }

        let vec_up_ = ifn.vecDiv(ifn.vecFromTo(x0_, x1_), num_y_ - 1);

        let vec_up_exp_ = ifn.vecSetLen(vec_up_, 100);

        for (let j_ of ifn.range(num_y_)) {

          let xyz_ = ifn.vecAdd(x0_, ifn.vecMult(vec_up_, j_));

          if (j_ == 0) {

            xyz_ = ifn.vecSub(xyz_, vec_up_exp_);
          } else {
            if (j_ == num_y_ - 1) {

              xyz_ = ifn.vecAdd(xyz_, vec_up_exp_);
            }
          }

          let posi_ = mfn.make.Position(xyz_);

          mfn.list.Add(posis_list_, posi_, 'to_end');
        }

        mfn.list.Add(posis_, posis_list_, 'to_end');
      }

      mfn.edit.Delete(quad_, 'delete_selected');

      return posis_;
    }


    async function exec_cluster1_node_gi0wr2wsff__createClusterParts_($p, off_grid_, off_grid_ex_, posis_) {

      let off_grid_e_ = mfn.query.Get('_e', off_grid_);

      let parts_ = [];

      let to_del_ = [];

      for (let i_ of ifn.range(ifn.len(posis_) - 1)) {

        let list0_ = posis_[pythonList(i_, posis_.length)];

        let list1_ = posis_[pythonList(i_ + 1, posis_.length)];

        for (let j_ of ifn.range(ifn.len(list0_) - 1)) {

          let posi0_ = list0_[pythonList(j_, list0_.length)];

          let posi1_ = list1_[pythonList(j_, list1_.length)];

          let posi2_ = list1_[pythonList(j_ + 1, list1_.length)];

          let posi3_ = list0_[pythonList(j_ + 1, list0_.length)];

          let part1_ = mfn.make.Polygon([posi0_, posi1_, posi2_, posi3_]);

          let part2_ = mfn.poly2d.Boolean(part1_, off_grid_ex_, 'intersect');

          await exec_cluster1_node_gi0wr2wsff__transferEdgeAttribsTouching_($p, off_grid_e_, mfn.query.Get('_e', part2_));
          if ($p.terminated) {
            return mfn.getModel();
          }

          mfn.list.Add(parts_, part2_, 'extend_end');

          mfn.list.Add(to_del_, part1_, 'to_end');
        }
      }

      mfn.edit.Delete(to_del_, 'delete_selected');

      return parts_;
    }


    async function exec_cluster1_node_gi0wr2wsff__trimTri_($p, off_grid_, tri_) {

      let tri2_ = mfn.poly2d.Boolean(tri_, off_grid_, 'intersect');

      mfn.edit.Delete(tri_, 'delete_selected');

      return tri2_;
    }


    async function exec_cluster1_node_gi0wr2wsff__vertAng_($p, vert_) {

      let edges_ = mfn.query.Get('_e', vert_);

      if (ifn.len(edges_) == 1) {

        return 0;
      }

      let vecs_ = mfn.calc.Vector(edges_);

      let vec0_ = ifn.vecRev(ifn.vecNorm(vecs_[pythonList(0, vecs_.length)]));

      let vec1_ = ifn.vecNorm(vecs_[pythonList(1, vecs_.length)]);

      return ifn.radToDeg(ifn.vecAng2(vec1_, vec0_, [0, 0, 1]));
    }


    async function exec_cluster1_node_gi0wr2wsff__angDot_($p, vert_) {

      let edges_ = mfn.query.Get('_e', vert_);

      let vecs_ = mfn.calc.Vector(edges_);

      let vec0_ = ifn.vecNorm(vecs_[pythonList(0, vecs_.length)]);

      let vec1_ = ifn.vecNorm(vecs_[pythonList(1, vecs_.length)]);

      return ifn.vecDot(vec0_, vec1_);
    }


    async function exec_cluster1_node_gi0wr2wsff__touchingEdge_($p, edges_, xyz_) {

      for (let edge_ of edges_) {

        let r_ = mfn.calc.Ray(edge_);

        let d_ = ifn.distance(xyz_, ifn.project(xyz_, r_, 0));

        if (d_ < 1) {

          return edge_;
        }
      }

      return null;
    }


    async function exec_cluster1_node_gi0wr2wsff__closestEdge_($p, edges_, xyz_) {

      let min_dist_ = Infinity;

      let close_edge_ = null;

      for (let edge_ of edges_) {

        let r_ = mfn.calc.Ray(edge_);

        let d_ = ifn.distance(xyz_, ifn.project(xyz_, r_, 0));

        if (d_ < min_dist_) {

          min_dist_ = d_;

          close_edge_ = edge_;
        }
      }

      return close_edge_;
    }


    async function exec_cluster1_node_gi0wr2wsff__transferEdgeAttribsTouching_($p, from_edges_, to_edges_) {

      for (let to_edge_ of to_edges_) {

        let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

        let from_edge_ = await exec_cluster1_node_gi0wr2wsff__touchingEdge_($p, from_edges_, cen_);
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


    async function exec_cluster1_node_gi0wr2wsff__transferEdgeAttribsClosest_($p, from_edges_, to_edges_) {

      for (let to_edge_ of to_edges_) {

        let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

        let from_edge_ = await exec_cluster1_node_gi0wr2wsff__closestEdge_($p, from_edges_, cen_);
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


    async function exec_cluster1_node_gi0wr2wsff__transferEdgeAttribsBtwTouchingParts_($p, parts_) {

      let edges_ = mfn.query.Filter(mfn.query.Get('_e', parts_), ['road', null], '==', null);

      for (let to_edge_ of edges_) {

        if (mfn.attrib.Get(to_edge_, 'road') == undefined) {

          let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

          let idx_ = ifn.listFind(edges_, to_edge_);

          let from_edges_ = ifn.listJoin(edges_.slice(0, idx_), edges_.slice(idx_ + 1));

          let from_edge_ = await exec_cluster1_node_gi0wr2wsff__touchingEdge_($p, from_edges_, cen_);
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


    async function exec_cluster1_node_gi0wr2wsff_isect_($p, ray_, rays_) {

      for (let a_ray_ of rays_) {

        let isect_ = ifn.intersect(ray_, a_ray_, 0);

        if (isect_ != null) {

          return isect_;
        }
      }

      return null;
    }

    async function exec_cluster1_node_gi0wr2wsff($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: subdiv into parts', '__null__')
      }


      let all_blocks_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', "block");

      mfn.edit.Delete(all_blocks_, 'delete_selected');

      let off_grids_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', "off_grid");

      for (let off_grid_ of off_grids_) {

        let block_type_ = mfn.attrib.Get(off_grid_, 'block_type');

        let type_ = mfn.attrib.Get(off_grid_, 'type');

        let parts_ = await exec_cluster1_node_gi0wr2wsff_subdivOfGrid_($p, off_grid_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (type_ != 'off_grid') {

          mfn.attrib.Set(parts_, `type`, type_);
        }

        mfn.attrib.Set(parts_, `block_type`, block_type_);
      }

      let on_grids_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '!=', "off_grid");

      for (let on_grid_ of on_grids_) {

        let block_type_ = mfn.attrib.Get(on_grid_, 'block_type');

        let type_ = mfn.attrib.Get(on_grid_, 'type');

        let parts_ = await exec_cluster1_node_gi0wr2wsff_subdivOfGrid_($p, on_grid_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (type_ != 'off_grid') {

          mfn.attrib.Set(parts_, `type`, type_);
        }

        mfn.attrib.Set(parts_, `block_type`, block_type_);
      }
    }


    async function exec_cluster1_node_xc8gbwkikvh_processBlock_($p, block_) {

      await exec_cluster1_node_xc8gbwkikvh__cleanPgonEdge_($p, block_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_node_xc8gbwkikvh__cleanPgonAng_($p, block_);
      if ($p.terminated) {
        return mfn.getModel();
      }
    }


    async function exec_cluster1_node_xc8gbwkikvh__angDot_($p, vert_) {

      let edges_ = mfn.query.Get('_e', vert_);

      let vecs_ = mfn.calc.Vector(edges_);

      let vec0_ = ifn.vecNorm(vecs_[pythonList(0, vecs_.length)]);

      let vec1_ = ifn.vecNorm(vecs_[pythonList(1, vecs_.length)]);

      return ifn.vecDot(vec0_, vec1_);
    }


    async function exec_cluster1_node_xc8gbwkikvh__cleanPgons_($p, pgons_) {

      for (let pgon_ of pgons_) {

        await exec_cluster1_node_xc8gbwkikvh__cleanPgonEdge_($p, pgon_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        await exec_cluster1_node_xc8gbwkikvh__cleanPgonAng_($p, pgon_);
        if ($p.terminated) {
          return mfn.getModel();
        }
      }

      return mfn.query.Get('pg', pgons_);
    }


    async function exec_cluster1_node_xc8gbwkikvh__cleanPgonEdge_($p, pgon_) {

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


    async function exec_cluster1_node_xc8gbwkikvh__cleanPgonAng_($p, pgon_) {

      let del_posis_ = [];

      for (let vert_ of mfn.query.Get('_v', pgon_)) {

        let dot_ = await exec_cluster1_node_xc8gbwkikvh__angDot_($p, vert_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (ifn.abs(dot_) > 0.9999) {

          mfn.list.Add(del_posis_, mfn.query.Get('ps', vert_), 'to_end');
        }
      }

      mfn.edit.Delete(del_posis_, 'delete_selected');
    }

    async function exec_cluster1_node_xc8gbwkikvh($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: cold blocks', '__null__')
      }


      let all_blocks_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', "block");

      let blocks0_ = mfn.query.Filter(mfn.query.Get('pg', all_blocks_), ['block_type', null], '==', "cold");

      let blocks_ = mfn.make.Clone(blocks0_);

      mfn.edit.Delete(blocks_, 'keep_selected');

      for (let block_ of blocks_) {

        await exec_cluster1_node_xc8gbwkikvh_processBlock_($p, block_);
        if ($p.terminated) {
          return mfn.getModel();
        }
      }
    }


    async function exec_cluster1_node_ytisso4wgrr_subdivBlock_($p, block_) {

      let block_id_ = mfn.attrib.Get(block_, 'block_id');

      let block_edges_ = mfn.query.Get('_e', block_);

      let block_verts_ = mfn.query.Get('_v', block_);

      let block_posis_ = mfn.query.Get('ps', block_verts_);

      let result_ = await exec_cluster1_node_ytisso4wgrr__getEdgesOnRoads_($p, block_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let road_edges_ = result_[pythonList(0, result_.length)];

      let other_edges_ = result_[pythonList(1, result_.length)];

      if (ifn.len(road_edges_) == 0) {

        return [];
      }

      let lengths_ = mfn.calc.Length(road_edges_);

      let edge_long_ = ifn.listSort(road_edges_, lengths_)[pythonList(-1, ifn.listSort(road_edges_, lengths_).length)];

      let vec_long_ = mfn.calc.Vector(edge_long_);

      vec_long_ = ifn.vecNorm(vec_long_);

      let vec_long_perp_ = [-vec_long_[pythonList(1, vec_long_.length)], vec_long_[pythonList(0, vec_long_.length)], 0];

      let cut0_ = null;

      let cut1_ = null;

      let posis0_ = null;

      let posis1_ = null;

      let subdivs_ = [];

      for (let i_ of ifn.range(1, ifn.len(road_edges_))) {

        cut0_ = cut1_;

        posis0_ = posis1_;

        let vert_ = mfn.query.Get('_v', road_edges_[pythonList(i_, road_edges_.length)])[pythonList(0, mfn.query.Get('_v', road_edges_[pythonList(i_, road_edges_.length)]).length)];

        let vert_edges_ = mfn.query.Get('_e', vert_);

        let road_idx0_ = ifn.listFind(block_edges_, vert_edges_[pythonList(0, vert_edges_.length)]);

        let road_idx1_ = ifn.listFind(block_edges_, vert_edges_[pythonList(1, vert_edges_.length)]);

        let posi_ = mfn.query.Get('ps', vert_)[pythonList(0, mfn.query.Get('ps', vert_).length)];

        let xyz_ = mfn.attrib.Get(posi_, 'xyz');

        let vecs_ = mfn.calc.Vector(vert_edges_);

        vecs_ = ifn.vecNorm(vecs_);

        let vec0_ = vecs_[pythonList(0, vecs_.length)];

        let vec1_ = vecs_[pythonList(1, vecs_.length)];

        let perp_vec0_ = [vec0_[pythonList(1, vec0_.length)], -vec0_[pythonList(0, vec0_.length)], 0];

        let perp_vec1_ = [vec1_[pythonList(1, vec1_.length)], -vec1_[pythonList(0, vec1_.length)], 0];

        let ang_ = await exec_cluster1_node_ytisso4wgrr__vertAng_($p, vert_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (ang_ > 250) {

          let xx_ = mfn.make.Extrude(vert_, 50, 1, 'quads');

          mfn.attrib.Set(xx_, `ang`, ang_);

          let xyz_a_ = null;

          let isect_a_ = [null, null];

          if (lengths_[pythonList(i_ - 1, lengths_.length)] > PART_LOC_D_) {

            let vec_a_ = ifn.vecRev(ifn.vecSetLen(perp_vec0_, 1000));

            xyz_a_ = ifn.vecAdd(xyz_, ifn.vecSetLen(vec0_, -(PLOT_LOC_W_ / 2)));

            let ray_a_ = ifn.rayMake(xyz_a_, vec_a_);

            isect_a_ = await exec_cluster1_node_ytisso4wgrr__isect_($p, ray_a_, other_edges_, subdivs_);
            if ($p.terminated) {
              return mfn.getModel();
            }

            if (isect_a_[pythonList(0, isect_a_.length)] != null && ifn.distance(xyz_a_, isect_a_[pythonList(0, isect_a_.length)]) > 300) {

              isect_a_ = [null, null];
            }
          }

          let xyz_b_ = null;

          let isect_b_ = [null, null];

          if (lengths_[pythonList(i_, lengths_.length)] > PART_LOC_D_) {

            let vec_b_ = ifn.vecRev(ifn.vecSetLen(perp_vec1_, 1000));

            xyz_b_ = ifn.vecAdd(xyz_, ifn.vecSetLen(vec1_, (PLOT_LOC_W_ / 2)));

            let ray_b_ = ifn.rayMake(xyz_b_, vec_b_);

            isect_b_ = await exec_cluster1_node_ytisso4wgrr__isect_($p, ray_b_, other_edges_, subdivs_);
            if ($p.terminated) {
              return mfn.getModel();
            }

            if (isect_b_[pythonList(0, isect_b_.length)] != null && ifn.distance(xyz_b_, isect_b_[pythonList(0, isect_b_.length)]) > 300) {

              isect_b_ = [null, null];
            }
          }

          let check_ = [isect_a_, isect_b_];

          if (isect_a_[pythonList(0, isect_a_.length)] != null && isect_b_[pythonList(0, isect_b_.length)] == null) {

            let vec_a_ = ifn.vecRev(ifn.vecSetLen(perp_vec0_, 1000));

            xyz_a_ = xyz_;

            let ray_a_ = ifn.rayMake(xyz_a_, vec_a_);

            isect_a_ = await exec_cluster1_node_ytisso4wgrr__isect_($p, ray_a_, other_edges_, subdivs_);
            if ($p.terminated) {
              return mfn.getModel();
            }
          }

          if (isect_a_[pythonList(0, isect_a_.length)] == null && isect_b_[pythonList(0, isect_b_.length)] != null) {

            let vec_b_ = ifn.vecRev(ifn.vecSetLen(perp_vec1_, 1000));

            xyz_b_ = xyz_;

            let ray_b_ = ifn.rayMake(xyz_b_, vec_b_);

            isect_b_ = await exec_cluster1_node_ytisso4wgrr__isect_($p, ray_b_, other_edges_, subdivs_);
            if ($p.terminated) {
              return mfn.getModel();
            }
          }

          if (isect_a_[pythonList(0, isect_a_.length)] == null && isect_b_[pythonList(0, isect_b_.length)] == null) {

            continue;
          } else {
            if (isect_a_[pythonList(0, isect_a_.length)] != null && isect_b_[pythonList(0, isect_b_.length)] == null) {

              posis1_ = mfn.make.Position([xyz_a_, isect_a_[pythonList(0, isect_a_.length)]]);

              let idx_a_ = ifn.listFind(block_edges_, isect_a_[pythonList(1, isect_a_.length)]);

              cut1_ = [road_idx0_, idx_a_];

              let subdiv_ = await exec_cluster1_node_ytisso4wgrr__makeSubdivs2_($p, block_posis_, cut0_, cut1_, posis0_, posis1_);
              if ($p.terminated) {
                return mfn.getModel();
              }

              mfn.attrib.Set(subdiv_, `type`, "block");

              mfn.list.Add(subdivs_, subdiv_, 'to_end');
            } else {
              if (isect_a_[pythonList(0, isect_a_.length)] == null && isect_b_[pythonList(0, isect_b_.length)] != null) {

                posis1_ = mfn.make.Position([xyz_b_, isect_b_[pythonList(0, isect_b_.length)]]);

                let idx_b_ = ifn.listFind(block_edges_, isect_b_[pythonList(1, isect_b_.length)]);

                cut1_ = [road_idx1_, idx_b_];

                let subdiv_ = await exec_cluster1_node_ytisso4wgrr__makeSubdivs2_($p, block_posis_, cut0_, cut1_, posis0_, posis1_);
                if ($p.terminated) {
                  return mfn.getModel();
                }

                mfn.attrib.Set(subdiv_, `type`, "block");

                mfn.list.Add(subdivs_, subdiv_, 'to_end');
              } else {

                let idx_a_ = ifn.listFind(block_edges_, isect_a_[pythonList(1, isect_a_.length)]);

                let dist_a_ = ifn.distance(isect_a_[pythonList(0, isect_a_.length)], xyz_);

                let idx_b_ = ifn.listFind(block_edges_, isect_b_[pythonList(1, isect_b_.length)]);

                let dist_b_ = ifn.distance(isect_b_[pythonList(0, isect_b_.length)], xyz_);

                posis1_ = mfn.make.Position([xyz_a_, isect_a_[pythonList(0, isect_a_.length)]]);

                cut1_ = [road_idx0_, idx_a_];

                let subdiv_ = await exec_cluster1_node_ytisso4wgrr__makeSubdivs2_($p, block_posis_, cut0_, cut1_, posis0_, posis1_);
                if ($p.terminated) {
                  return mfn.getModel();
                }

                mfn.attrib.Set(subdiv_, `type`, "block");

                mfn.list.Add(subdivs_, subdiv_, 'to_end');

                posis0_ = posis1_;

                cut0_ = cut1_;

                posis1_ = mfn.make.Position([xyz_b_, isect_b_[pythonList(0, isect_b_.length)]]);

                cut1_ = [road_idx1_, idx_b_];

                subdiv_ = await exec_cluster1_node_ytisso4wgrr__makeSubdivs2_($p, block_posis_, cut0_, cut1_, posis0_, posis1_);
                if ($p.terminated) {
                  return mfn.getModel();
                }

                mfn.attrib.Set(subdiv_, `type`, "block_corner");

                mfn.visualize.Color(subdiv_, [0.7, 0, 0.7]);

                mfn.list.Add(subdivs_, subdiv_, 'to_end');
              }
            }
          }
        } else {

          continue;
        }
      }

      posis0_ = posis1_;

      cut0_ = cut1_;

      cut1_ = null;

      let subdiv_ = await exec_cluster1_node_ytisso4wgrr__makeSubdivs2_($p, block_posis_, cut0_, cut1_, posis0_, posis1_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      mfn.attrib.Set(subdiv_, `type`, "block");

      mfn.list.Add(subdivs_, subdiv_, 'to_end');

      await exec_cluster1_node_ytisso4wgrr__transferEdgeAttribs_($p, block_edges_, mfn.query.Get('_e', subdivs_));
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_node_ytisso4wgrr__transferEdgeAttribsBtwTouchingParts_($p, subdivs_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_node_ytisso4wgrr__copyAttribs_($p, block_, subdivs_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      return subdivs_;
    }


    async function exec_cluster1_node_ytisso4wgrr__transferEdgeAttribs_($p, from_edges_, to_edges_) {

      for (let to_edge_ of to_edges_) {

        let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

        let from_edge_ = await exec_cluster1_node_ytisso4wgrr__touchingEdge_($p, from_edges_, cen_);
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


    async function exec_cluster1_node_ytisso4wgrr__transferEdgeAttribsBtwTouchingParts_($p, parts_) {

      let edges_ = mfn.query.Filter(mfn.query.Get('_e', parts_), ['road', null], '==', null);

      for (let to_edge_ of edges_) {

        if (mfn.attrib.Get(to_edge_, 'road') == undefined) {

          let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

          let idx_ = ifn.listFind(edges_, to_edge_);

          let from_edges_ = ifn.listJoin(edges_.slice(0, idx_), edges_.slice(idx_ + 1));

          let from_edge_ = await exec_cluster1_node_ytisso4wgrr__touchingEdge_($p, from_edges_, cen_);
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


    async function exec_cluster1_node_ytisso4wgrr__copyAttribs_($p, pgon_from_, pgons_to_) {

      let block_id_ = mfn.attrib.Get(pgon_from_, 'block_id');

      let block_type_ = mfn.attrib.Get(pgon_from_, 'block_type');

      let site_ = mfn.attrib.Get(pgon_from_, 'site');

      for (let pgon_to_ of pgons_to_) {

        mfn.attrib.Set(pgon_to_, `block_id`, block_id_);

        mfn.attrib.Set(pgon_to_, `block_type`, block_type_);

        mfn.attrib.Set(pgon_to_, `site`, site_);
      }
    }


    async function exec_cluster1_node_ytisso4wgrr__touchingEdge_($p, edges_, xyz_) {

      for (let edge_ of edges_) {

        let r_ = mfn.calc.Ray(edge_);

        let d_ = ifn.distance(xyz_, ifn.project(xyz_, r_, 0));

        if (d_ < 1) {

          return edge_;
        }
      }

      return null;
    }


    async function exec_cluster1_node_ytisso4wgrr__getEdgesOnRoads_($p, block_) {

      let road_edges_lists_ = [[]];

      let other_edges_ = [];

      for (let edge_ of mfn.query.Get('_e', block_)) {

        let edge_road_ = mfn.attrib.Get(edge_, 'road');

        if (edge_road_ == "road_art" || edge_road_ == "road_loc" || edge_road_ == "road_sec") {

          mfn.list.Add(road_edges_lists_[pythonList(-1, road_edges_lists_.length)], edge_, 'to_end');
        } else {

          mfn.list.Add(other_edges_, edge_, 'to_end');

          if (ifn.len(road_edges_lists_[pythonList(-1, road_edges_lists_.length)]) != 0) {

            mfn.list.Add(road_edges_lists_, [], 'to_end');
          }
        }
      }

      let road_edges_ = [];

      if (ifn.len(road_edges_lists_) == 0) {

        road_edges_ = [];
      } else {
        if (ifn.len(road_edges_lists_) == 1) {

          road_edges_ = road_edges_lists_[pythonList(0, road_edges_lists_.length)];
        } else {

          road_edges_ = ifn.listJoin(road_edges_lists_[pythonList(1, road_edges_lists_.length)], road_edges_lists_[pythonList(0, road_edges_lists_.length)]);
        }
      }

      return [road_edges_, other_edges_];
    }


    async function exec_cluster1_node_ytisso4wgrr__isect_($p, ray_, edges_, subdivs_) {

      let rays_ = mfn.calc.Ray(edges_);

      for (let edge_ of mfn.query.Get('_e', subdivs_)) {

        let edge_ray_ = mfn.calc.Ray(edge_);

        let isect_ = ifn.intersect(ray_, edge_ray_, 0);

        if (isect_ != null) {

          return [null, null];
        }
      }

      let isect_min_ = null;

      let dist_min_ = Infinity;

      let edge_min_ = null;

      for (let edge_ of edges_) {

        let edge_ray_ = mfn.calc.Ray(edge_);

        let isect_ = ifn.intersect(ray_, edge_ray_, 0);

        if (isect_ != null) {

          let dist_ = ifn.distance(isect_, ray_[pythonList(0, ray_.length)]);

          if (dist_ < dist_min_) {

            isect_min_ = isect_;

            dist_min_ = dist_;

            edge_min_ = edge_;
          }
        }
      }

      return [isect_min_, edge_min_];
    }


    async function exec_cluster1_node_ytisso4wgrr__vertAng_($p, vert_) {

      let edges_ = mfn.query.Get('_e', vert_);

      if (ifn.len(edges_) == 1) {

        return 0;
      }

      let vecs_ = mfn.calc.Vector(edges_);

      let vec0_ = ifn.vecRev(ifn.vecNorm(vecs_[pythonList(0, vecs_.length)]));

      let vec1_ = ifn.vecNorm(vecs_[pythonList(1, vecs_.length)]);

      return ifn.radToDeg(ifn.vecAng2(vec1_, vec0_, [0, 0, 1]));
    }


    async function exec_cluster1_node_ytisso4wgrr__makeSubdivs2_($p, block_posis_, cut0_, cut1_, posis0_, posis1_) {

      let ring_ = [];

      if (cut0_ == null && cut1_ == null) {

        let pgon_ = mfn.make.Polygon(block_posis_);

        mfn.attrib.Set(pgon_, `type`, "off_grid");

        return pgon_;
      } else {
        if (cut0_ == null) {

          ring_ = ifn.listCopy(posis1_);

          let posis_list_ = await exec_cluster1_node_ytisso4wgrr__getPosisFromRing_($p, block_posis_, [cut1_[pythonList(1, cut1_.length)], cut1_[pythonList(0, cut1_.length)]]);
          if ($p.terminated) {
            return mfn.getModel();
          }

          mfn.list.Add(ring_, posis_list_, 'extend_end');
        } else {
          if (cut1_ == null) {

            ring_ = ifn.listRev(posis0_);

            let posis_list_ = await exec_cluster1_node_ytisso4wgrr__getPosisFromRing_($p, block_posis_, [cut0_[pythonList(0, cut0_.length)], cut0_[pythonList(1, cut0_.length)]]);
            if ($p.terminated) {
              return mfn.getModel();
            }

            mfn.list.Add(ring_, posis_list_, 'extend_end');
          } else {

            ring_ = ifn.listCopy(posis1_);

            let posis_list_ = await exec_cluster1_node_ytisso4wgrr__getPosisFromRing_($p, block_posis_, [cut1_[pythonList(1, cut1_.length)], cut0_[pythonList(1, cut0_.length)]]);
            if ($p.terminated) {
              return mfn.getModel();
            }

            mfn.list.Add(ring_, posis_list_, 'extend_end');

            mfn.list.Add(ring_, ifn.listRev(posis0_), 'extend_end');

            posis_list_ = await exec_cluster1_node_ytisso4wgrr__getPosisFromRing_($p, block_posis_, [cut0_[pythonList(0, cut0_.length)], cut1_[pythonList(0, cut1_.length)]]);
            if ($p.terminated) {
              return mfn.getModel();
            }

            mfn.list.Add(ring_, posis_list_, 'extend_end');
          }
        }
      }

      if (ifn.len(ring_) < 3) {

        return null;
      }

      let pgon_ = mfn.make.Polygon(ring_);

      return pgon_;
    }


    async function exec_cluster1_node_ytisso4wgrr__getPosisFromRing_($p, posis_, idxs_) {

      let num_posis_ = ifn.len(posis_);

      let idx0_ = (idxs_[pythonList(0, idxs_.length)] + 1) % num_posis_;

      let idx1_ = (idxs_[pythonList(1, idxs_.length)] + 1) % num_posis_;

      if (idx0_ == idx1_) {

        return [];
      }

      let ring_ = [];

      if (idx0_ < idx1_) {

        mfn.list.Add(ring_, posis_.slice(idx0_, idx1_), 'extend_end');
      } else {

        mfn.list.Add(ring_, posis_.slice(idx0_), 'extend_end');

        mfn.list.Add(ring_, posis_.slice(0, idx1_), 'extend_end');
      }

      return ring_;
    }

    async function exec_cluster1_node_ytisso4wgrr($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: sub div cold blocks at concave corners', '__null__')
      }


      let blocks_ = mfn.query.Get('pg', null);

      for (let block_ of blocks_) {

        let subdivs_ = await exec_cluster1_node_ytisso4wgrr_subdivBlock_($p, block_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        let count_ = 0;

        for (let subdiv_ of subdivs_) {

          mfn.attrib.Set(subdiv_, `block_id`, mfn.attrib.Get(subdiv_, 'block_id') + "_" + count_);

          count_ = count_ + 1;
        }
      }

      mfn.edit.Delete([blocks_], 'delete_selected');
    }


    async function exec_cluster1_node_joqspe73iin($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: get blocks', '__null__')
      }


      let blocks_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', "block");

      mfn.edit.Delete(blocks_, 'keep_selected');
    }


    async function exec_cluster1_node_bzwyrx9mkk6_processBlock_($p, block_) {

      let roads_art_ = await exec_cluster1_node_bzwyrx9mkk6__getSitePlines_($p, block_, "road_art");
      if ($p.terminated) {
        return mfn.getModel();
      }

      let roads_sec_ = await exec_cluster1_node_bzwyrx9mkk6__getSitePlines_($p, block_, "road_sec");
      if ($p.terminated) {
        return mfn.getModel();
      }

      let roads_loc_ = await exec_cluster1_node_bzwyrx9mkk6__getSitePlines_($p, block_, "road_loc");
      if ($p.terminated) {
        return mfn.getModel();
      }

      let parts_ = await exec_cluster1_node_bzwyrx9mkk6__createOnGridStrips_($p, block_, roads_art_, roads_sec_, roads_loc_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      parts_ = await exec_cluster1_node_bzwyrx9mkk6__cleanPgons_($p, parts_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      return parts_;
    }


    async function exec_cluster1_node_bzwyrx9mkk6__createOnGridStrips_($p, block_, roads_art_, roads_sec_, roads_loc_) {

      let off_art_ = mfn.poly2d.OffsetMitre(roads_art_, PART_ART_D_, 100, 'square_end');

      let off_sec_ = mfn.poly2d.OffsetMitre(roads_sec_, PART_SEC_D_, 100, 'square_end');

      let off_loc_ = mfn.poly2d.OffsetMitre(roads_loc_, PART_LOC_D_, 100, 'square_end');

      let off_grids_ = mfn.poly2d.Boolean(block_, [off_art_, off_sec_, off_loc_], 'difference');

      let on_arts1_ = mfn.poly2d.Boolean(off_art_, block_, 'intersect');

      let on_secs1_ = mfn.poly2d.Boolean(off_sec_, block_, 'intersect');

      let on_locs1_ = mfn.poly2d.Boolean(off_loc_, block_, 'intersect');

      let on_arts2_ = mfn.poly2d.Boolean(on_arts1_, [off_sec_, off_loc_], 'difference');

      let on_secs2_ = mfn.poly2d.Boolean(on_secs1_, [off_art_, off_loc_], 'difference');

      let on_locs2_ = mfn.poly2d.Boolean(on_locs1_, [off_art_, off_sec_], 'difference');

      let on_arts3_ = mfn.poly2d.Union(on_arts2_);

      let on_secs3_ = mfn.poly2d.Union(on_secs2_);

      let on_locs3_ = mfn.poly2d.Union(on_locs2_);

      let art_sec_ = mfn.poly2d.Boolean(off_art_, off_sec_, 'intersect');

      let art_loc_ = mfn.poly2d.Boolean(off_art_, off_loc_, 'intersect');

      let sec_loc_ = mfn.poly2d.Boolean(off_sec_, off_loc_, 'intersect');

      let art_sec1_ = mfn.poly2d.Boolean(art_sec_, block_, 'intersect');

      let art_loc1_ = mfn.poly2d.Boolean(art_loc_, block_, 'intersect');

      let sec_loc1_ = mfn.poly2d.Boolean(sec_loc_, block_, 'intersect');

      let new_parts_ = ifn.listFlat([art_sec1_, art_loc1_, sec_loc1_, on_arts3_, on_secs3_, on_locs3_, off_grids_]);

      mfn.attrib.Set(art_sec1_, `type`, 'art_sec');

      mfn.attrib.Set(art_loc1_, `type`, 'art_loc');

      mfn.attrib.Set(sec_loc1_, `type`, 'sec_loc');

      mfn.attrib.Set(on_arts3_, `type`, 'art');

      mfn.attrib.Set(on_secs3_, `type`, 'sec');

      mfn.attrib.Set(on_locs3_, `type`, 'loc');

      mfn.attrib.Set(off_grids_, `type`, 'off_grid');

      mfn.attrib.Set(new_parts_, `class`, "part");

      await exec_cluster1_node_bzwyrx9mkk6__transferEdgeAttribsBtwTouchingParts_($p, new_parts_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_node_bzwyrx9mkk6__transferEdgeAttribsTouching_($p, mfn.query.Get('_e', block_), mfn.query.Filter(mfn.query.Get('_e', new_parts_), ['road', null], '==', null));
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_node_bzwyrx9mkk6__copyAttribs_($p, block_, new_parts_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      mfn.edit.Delete([art_sec_, art_loc_, sec_loc_], 'delete_selected');

      mfn.edit.Delete([off_art_, off_sec_, off_loc_, on_arts1_, on_arts2_, on_secs1_, on_secs2_, on_locs1_, on_locs2_, block_], 'delete_selected');

      mfn.edit.Delete([roads_art_, roads_sec_, roads_loc_], 'delete_selected');

      return new_parts_;
    }


    async function exec_cluster1_node_bzwyrx9mkk6__angDot_($p, vert_) {

      let edges_ = mfn.query.Get('_e', vert_);

      let vecs_ = mfn.calc.Vector(edges_);

      let vec0_ = ifn.vecNorm(vecs_[pythonList(0, vecs_.length)]);

      let vec1_ = ifn.vecNorm(vecs_[pythonList(1, vecs_.length)]);

      return ifn.vecDot(vec0_, vec1_);
    }


    async function exec_cluster1_node_bzwyrx9mkk6__cleanPgons_($p, pgons_) {

      for (let pgon_ of pgons_) {

        await exec_cluster1_node_bzwyrx9mkk6__cleanPgonEdge_($p, pgon_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        await exec_cluster1_node_bzwyrx9mkk6__cleanPgonAng_($p, pgon_);
        if ($p.terminated) {
          return mfn.getModel();
        }
      }

      return mfn.query.Get('pg', pgons_);
    }


    async function exec_cluster1_node_bzwyrx9mkk6__cleanPgonEdge_($p, pgon_) {

      let del_posis_ = [];

      for (let edge_ of mfn.query.Get('_e', pgon_)) {

        let length_ = mfn.calc.Length(edge_);

        if (length_ < 0.0001) {

          let posis_ = mfn.query.Get('ps', edge_);

          mfn.list.Add(del_posis_, posis_[pythonList(0, posis_.length)], 'to_end');
        }
      }

      mfn.edit.Delete(del_posis_, 'delete_selected');
    }


    async function exec_cluster1_node_bzwyrx9mkk6__cleanPgonAng_($p, pgon_) {

      let del_posis_ = [];

      for (let vert_ of mfn.query.Get('_v', pgon_)) {

        let dot_ = await exec_cluster1_node_bzwyrx9mkk6__angDot_($p, vert_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (ifn.abs(dot_) > 0.9999) {

          mfn.list.Add(del_posis_, mfn.query.Get('ps', vert_), 'to_end');
        }
      }

      mfn.edit.Delete(del_posis_, 'delete_selected');
    }


    async function exec_cluster1_node_bzwyrx9mkk6__getSitePlines_($p, site_, road_descr_) {

      let posis_ = [];

      for (let edge_ of mfn.query.Get('_e', site_)) {

        if (mfn.attrib.Get(edge_, 'road') == road_descr_) {

          let start_end_ = mfn.query.Get('ps', edge_);

          if (ifn.len(posis_) == 0 || start_end_[pythonList(0, start_end_.length)] != posis_[pythonList(-1, posis_.length)][pythonList(-1, posis_[pythonList(-1, posis_.length)].length)]) {

            mfn.list.Add(posis_, start_end_, 'to_end');
          } else {

            mfn.list.Add(posis_[pythonList(-1, posis_.length)], start_end_[pythonList(1, start_end_.length)], 'to_end');
          }
        }
      }

      if (ifn.len(posis_) == 0) {

        return [];
      }

      if (ifn.len(posis_) > 1 && posis_[pythonList(-1, posis_.length)][pythonList(-1, posis_[pythonList(-1, posis_.length)].length)] == posis_[pythonList(0, posis_.length)][pythonList(0, posis_[pythonList(0, posis_.length)].length)]) {

        let first_list_ = ifn.listJoin(posis_[pythonList(-1, posis_.length)], posis_[pythonList(0, posis_.length)].slice(1));

        posis_[pythonList(0, posis_.length)] = first_list_;

        posis_ = posis_.slice(0, -1);
      }

      let site_plines_ = mfn.make.Polyline(posis_, 'open');

      return site_plines_;
    }


    async function exec_cluster1_node_bzwyrx9mkk6__transferEdgeAttribsTouching_($p, from_edges_, to_edges_) {

      for (let to_edge_ of to_edges_) {

        let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

        let from_edge_ = await exec_cluster1_node_bzwyrx9mkk6__touchingEdge_($p, from_edges_, cen_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (from_edge_ != null) {

          mfn.attrib.Set(to_edge_, `road`, mfn.attrib.Get(from_edge_, 'road'));
        }
      }
    }


    async function exec_cluster1_node_bzwyrx9mkk6__transferEdgeAttribsBtwTouchingParts_($p, parts_) {

      let edges_ = mfn.query.Filter(mfn.query.Get('_e', parts_), ['road', null], '==', null);

      for (let to_edge_ of edges_) {

        if (mfn.attrib.Get(to_edge_, 'road') == undefined) {

          let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

          let idx_ = ifn.listFind(edges_, to_edge_);

          let from_edges_ = ifn.listJoin(edges_.slice(0, idx_), edges_.slice(idx_ + 1));

          let from_edge_ = await exec_cluster1_node_bzwyrx9mkk6__touchingEdge_($p, from_edges_, cen_);
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


    async function exec_cluster1_node_bzwyrx9mkk6__copyAttribs_($p, pgon_from_, pgons_to_) {

      let block_id_ = mfn.attrib.Get(pgon_from_, 'block_id');

      let block_type_ = mfn.attrib.Get(pgon_from_, 'block_type');

      let site_ = mfn.attrib.Get(pgon_from_, 'site');

      for (let pgon_to_ of pgons_to_) {

        mfn.attrib.Set(pgon_to_, `block_id`, block_id_);

        mfn.attrib.Set(pgon_to_, `block_type`, block_type_);

        mfn.attrib.Set(pgon_to_, `site`, site_);
      }
    }


    async function exec_cluster1_node_bzwyrx9mkk6__touchingEdge_($p, edges_, xyz_) {

      for (let edge_ of edges_) {

        let r_ = mfn.calc.Ray(edge_);

        let d_ = ifn.distance(xyz_, ifn.project(xyz_, r_, 0));

        if (d_ < 0.01) {

          return edge_;
        }
      }

      return null;
    }

    async function exec_cluster1_node_bzwyrx9mkk6($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: add on grid strips art sec loc', '__null__')
      }


      let blocks_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', "block");

      mfn.edit.Delete(blocks_, 'keep_selected');

      for (let block_ of blocks_) {

        let result_ = await exec_cluster1_node_bzwyrx9mkk6_processBlock_($p, block_);
        if ($p.terminated) {
          return mfn.getModel();
        }
      }
    }


    async function exec_cluster1_node_q0cbraoo5qe_subdivParts_($p, part_) {

      let part_type_ = mfn.attrib.Get(part_, 'type');

      let part_depth_ = PART_SEC_D_;

      let plot_width_ = PLOT_SEC_W_;

      if (part_type_ == "loc") {

        part_depth_ = PART_LOC_D_;

        plot_width_ = PLOT_LOC_W_;
      }

      let part_edges_ = mfn.query.Get('_e', part_);

      let road_loc_e_ = mfn.query.Filter(part_edges_, ['road', null], '==', "road_loc");

      let road_sec_e_ = mfn.query.Filter(part_edges_, ['road', null], '==', "road_sec");

      let road_art_e_ = mfn.query.Filter(part_edges_, ['road', null], '==', "road_art");

      let og_e_ = mfn.query.Filter(part_edges_, ['road', null], '==', "off_grid");

      let road_edges_ = ifn.listFlat([road_loc_e_, road_sec_e_, road_art_e_]);

      let ang_threshold_ = 200;

      let og_verts_ = mfn.query.Get('_v', og_e_);

      let convex_ = [];

      for (let vert_ of og_verts_) {

        let ang_ = await exec_cluster1_node_q0cbraoo5qe__vertAng_($p, vert_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (ang_ > ang_threshold_) {

          let len_edges_ = mfn.calc.Length(mfn.query.Get('_e', vert_));

          if (len_edges_[pythonList(0, len_edges_.length)] > (PLOT_LOC_W_ / 2) && len_edges_[pythonList(1, len_edges_.length)] > (PLOT_LOC_W_ / 2)) {

            let corner_ = await exec_cluster1_node_q0cbraoo5qe__makeConvexCorner_($p, vert_, part_depth_);
            if ($p.terminated) {
              return mfn.getModel();
            }

            mfn.list.Add(convex_, corner_, 'to_end');
          }
        }
      }

      let corners_un_ = mfn.poly2d.Union(convex_);

      let corners_trim_ = mfn.poly2d.Boolean(corners_un_, part_, 'intersect');

      let strips_ = mfn.poly2d.Boolean(part_, corners_un_, 'difference');

      strips_ = await exec_cluster1_node_q0cbraoo5qe__fixPinchedPgons_($p, strips_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_node_q0cbraoo5qe__attribs_($p, part_, corners_trim_, strips_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      mfn.edit.Delete([part_, convex_, corners_un_], 'delete_selected');

      return [strips_, corners_trim_];
    }


    async function exec_cluster1_node_q0cbraoo5qe__attribs_($p, part_, corners_, strips_) {

      let part_edges_ = mfn.query.Get('_e', part_);

      let part_type_ = mfn.attrib.Get(part_, 'type');

      await exec_cluster1_node_q0cbraoo5qe__transferEdgeAttribsTouching_($p, part_edges_, mfn.query.Get('_e', corners_));
      if ($p.terminated) {
        return mfn.getModel();
      }

      mfn.attrib.Set(corners_, `type`, part_type_ + '_' + part_type_);

      await exec_cluster1_node_q0cbraoo5qe__transferEdgeAttribsTouching_($p, part_edges_, mfn.query.Get('_e', strips_));
      if ($p.terminated) {
        return mfn.getModel();
      }

      mfn.attrib.Set(strips_, `type`, part_type_);

      let all_parts_ = ifn.listFlat([strips_, corners_]);

      mfn.attrib.Set(all_parts_, `class`, "part");

      await exec_cluster1_node_q0cbraoo5qe__copyAttribs_($p, part_, all_parts_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_node_q0cbraoo5qe__transferEdgeAttribsBtwTouchingParts_($p, all_parts_);
      if ($p.terminated) {
        return mfn.getModel();
      }
    }


    async function exec_cluster1_node_q0cbraoo5qe__makeConvexCorner_($p, vert_, part_depth_) {

      let vert_xyz_ = mfn.attrib.Get(mfn.query.Get('ps', vert_)[pythonList(0, mfn.query.Get('ps', vert_).length)], 'xyz');

      let vert_edges_ = mfn.query.Get('_e', vert_);

      let vec0_ = mfn.calc.Vector(vert_edges_[pythonList(0, vert_edges_.length)]);

      let vec1_ = mfn.calc.Vector(vert_edges_[pythonList(1, vert_edges_.length)]);

      let vec0_perp_ = ifn.vecSetLen([-vec0_[pythonList(1, vec0_.length)], vec0_[pythonList(0, vec0_.length)], 0], part_depth_ + 1);

      let vec1_perp_ = ifn.vecSetLen([-vec1_[pythonList(1, vec1_.length)], vec1_[pythonList(0, vec1_.length)], 0], part_depth_ + 1);

      let xyz0_ = ifn.vecAdd(vert_xyz_, vec0_perp_);

      let xyz1_ = ifn.vecAdd(vert_xyz_, vec1_perp_);

      let ray0_ = ifn.rayMake(xyz0_, vec0_);

      let ray1_ = ifn.rayMake(xyz1_, vec1_);

      let xyz_cor_ = ifn.intersect(ray0_, ray1_, 2);

      let posis_ = mfn.make.Position([xyz1_, xyz_cor_, xyz0_]);

      posis_ = ifn.listJoin(mfn.query.Get('ps', vert_), posis_);

      let pgon_ = mfn.make.Polygon(posis_);

      return pgon_;
    }


    async function exec_cluster1_node_q0cbraoo5qe__makeConcaveCorner_($p, vert_, part_depth_, plot_width_) {

      let vert_xyz_ = mfn.attrib.Get(mfn.query.Get('ps', vert_)[pythonList(0, mfn.query.Get('ps', vert_).length)], 'xyz');

      let vert_edges_ = mfn.query.Get('_e', vert_);

      let vec0_ = mfn.calc.Vector(vert_edges_[pythonList(0, vert_edges_.length)]);

      let vec1_ = mfn.calc.Vector(vert_edges_[pythonList(1, vert_edges_.length)]);

      let corner_front_ = plot_width_ / 2;

      let xyz_cor_off_ = ifn.vecSum(vert_xyz_, ifn.vecSetLen(vec0_, -corner_front_), ifn.vecSetLen(vec1_, corner_front_));

      let vec0_perp_ = ifn.vecSetLen([-vec0_[pythonList(1, vec0_.length)], vec0_[pythonList(0, vec0_.length)], 0], part_depth_ + corner_front_ + 1);

      let vec1_perp_ = ifn.vecSetLen([-vec1_[pythonList(1, vec1_.length)], vec1_[pythonList(0, vec1_.length)], 0], part_depth_ + corner_front_ + 1);

      let xyz0_ = ifn.vecAdd(xyz_cor_off_, vec0_perp_);

      let xyz1_ = ifn.vecAdd(xyz_cor_off_, vec1_perp_);

      let ray0_ = ifn.rayMake(xyz0_, vec0_);

      let ray1_ = ifn.rayMake(xyz1_, vec1_);

      let xyz_cor_ = ifn.intersect(ray0_, ray1_, 2);

      let posis_ = mfn.make.Position([xyz_cor_off_, xyz1_, xyz_cor_, xyz0_]);

      let pgon_ = mfn.make.Polygon(posis_);

      return pgon_;
    }


    async function exec_cluster1_node_q0cbraoo5qe__makeZigZagCorner_($p, edge_, part_depth_, plot_width_) {

      let verts_ = mfn.query.Get('_v', edge_);

      let vec0_ = mfn.calc.Vector(edge_);

      let length_ = ifn.vecLen(vec0_);

      vec0_ = ifn.vecNorm(vec0_);

      let vec_next_ = mfn.calc.Vector(mfn.query.Get('_e', verts_[pythonList(1, verts_.length)])[pythonList(1, mfn.query.Get('_e', verts_[pythonList(1, verts_.length)]).length)]);

      let vec_next_perp_ = [-vec_next_[pythonList(1, vec_next_.length)], vec_next_[pythonList(0, vec_next_.length)], 0];

      let vec_prev_ = mfn.calc.Vector(mfn.query.Get('_e', verts_[pythonList(0, verts_.length)])[pythonList(0, mfn.query.Get('_e', verts_[pythonList(0, verts_.length)]).length)]);

      let vec_prev_perp_ = [-vec_prev_[pythonList(1, vec_prev_.length)], vec_prev_[pythonList(0, vec_prev_.length)], 0];

      let posi_ = mfn.query.Get('ps', verts_[pythonList(0, verts_.length)])[pythonList(0, mfn.query.Get('ps', verts_[pythonList(0, verts_.length)]).length)];

      let xyz0_ = mfn.attrib.Get(posi_, 'xyz');

      let xyz1_ = ifn.vecAdd(xyz0_, ifn.vecSetLen(vec0_, length_ + 1));

      let xyz2_ = ifn.vecAdd(xyz1_, ifn.vecSetLen(vec_next_, part_depth_ + plot_width_));

      let xyz3_ = ifn.vecAdd(xyz2_, ifn.vecSetLen(vec_next_perp_, length_ + part_depth_ + 10));

      let xyz4_ = ifn.vecAdd(xyz0_, ifn.vecSetLen(vec_prev_perp_, part_depth_ + 10));

      let posis_ = mfn.make.Position([xyz1_, xyz2_, xyz3_, xyz4_]);

      mfn.list.Add(posis_, posi_, 'extend_end');

      let pgon_ = mfn.make.Polygon(posis_);

      return pgon_;
    }


    async function exec_cluster1_node_q0cbraoo5qe__fixPinchedPgons_($p, pgons_) {

      let all_new_pgons_ = [];

      for (let pgon_ of pgons_) {

        let verts_ = mfn.query.Get('_v', pgon_);

        let vert_posis_ = [];

        let pinch_posis_ = [];

        for (let vert_ of verts_) {

          let posi_ = mfn.query.Get('ps', vert_)[pythonList(0, mfn.query.Get('ps', vert_).length)];

          if (ifn.listHas(vert_posis_, posi_)) {

            mfn.list.Add(pinch_posis_, posi_, 'to_end');
          }

          mfn.list.Add(vert_posis_, posi_, 'to_end');
        }

        let new_posis_ = [[]];

        let idx_ = -1;

        let used_ = [];

        for (let posi_ of vert_posis_) {

          if (ifn.listHas(pinch_posis_, posi_)) {

            if (ifn.listHas(used_, posi_)) {

              idx_ = idx_ - 1;
            } else {

              mfn.list.Add(used_, posi_, 'to_end');

              mfn.list.Add(new_posis_, [], 'to_end');
            }
          }

          mfn.list.Add(new_posis_[pythonList(idx_, new_posis_.length)], posi_, 'to_end');
        }

        for (let a_new_posis_ of new_posis_) {

          if (ifn.len(a_new_posis_) > 2) {

            let new_pgon_ = mfn.make.Polygon(a_new_posis_);

            mfn.list.Add(all_new_pgons_, new_pgon_, 'extend_end');
          }
        }
      }

      mfn.edit.Delete(pgons_, 'delete_selected');

      return all_new_pgons_;
    }


    async function exec_cluster1_node_q0cbraoo5qe__vertAng_($p, vert_) {

      let edges_ = mfn.query.Get('_e', vert_);

      if (ifn.len(edges_) == 1) {

        return 0;
      }

      let vecs_ = mfn.calc.Vector(edges_);

      let vec0_ = ifn.vecRev(ifn.vecNorm(vecs_[pythonList(0, vecs_.length)]));

      let vec1_ = ifn.vecNorm(vecs_[pythonList(1, vecs_.length)]);

      return ifn.radToDeg(ifn.vecAng2(vec1_, vec0_, [0, 0, 1]));
    }


    async function exec_cluster1_node_q0cbraoo5qe__transferEdgeAttribsTouching_($p, from_edges_, to_edges_) {

      for (let to_edge_ of to_edges_) {

        let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

        let from_edge_ = await exec_cluster1_node_q0cbraoo5qe__touchingEdge_($p, from_edges_, cen_);
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


    async function exec_cluster1_node_q0cbraoo5qe__transferEdgeAttribsTouchingPart_($p, from_edges_, to_edges_) {

      for (let to_edge_ of to_edges_) {

        if (mfn.attrib.Get(to_edge_, 'road') == undefined) {

          let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

          let from_edge_ = await exec_cluster1_node_q0cbraoo5qe__touchingEdge_($p, from_edges_, cen_);
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


    async function exec_cluster1_node_q0cbraoo5qe__transferEdgeAttribsBtwTouchingParts_($p, parts_) {

      let edges_ = mfn.query.Get('_e', parts_);

      for (let to_edge_ of edges_) {

        if (mfn.attrib.Get(to_edge_, 'road') == undefined) {

          let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

          let idx_ = ifn.listFind(edges_, to_edge_);

          let from_edges_ = ifn.listJoin(edges_.slice(0, idx_), edges_.slice(idx_ + 1));

          let from_edge_ = await exec_cluster1_node_q0cbraoo5qe__touchingEdge_($p, from_edges_, cen_);
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


    async function exec_cluster1_node_q0cbraoo5qe__copyAttribs_($p, pgon_from_, pgons_to_) {

      let block_id_ = mfn.attrib.Get(pgon_from_, 'block_id');

      let block_type_ = mfn.attrib.Get(pgon_from_, 'block_type');

      let site_ = mfn.attrib.Get(pgon_from_, 'site');

      for (let pgon_to_ of pgons_to_) {

        mfn.attrib.Set(pgon_to_, `block_id`, block_id_);

        mfn.attrib.Set(pgon_to_, `block_type`, block_type_);

        mfn.attrib.Set(pgon_to_, `site`, site_);
      }
    }


    async function exec_cluster1_node_q0cbraoo5qe__touchingEdge_($p, edges_, xyz_) {

      for (let edge_ of edges_) {

        let r_ = mfn.calc.Ray(edge_);

        let d_ = ifn.distance(xyz_, ifn.project(xyz_, r_, 0));

        if (d_ < 0.1) {

          return edge_;
        }
      }

      return null;
    }

    async function exec_cluster1_node_q0cbraoo5qe($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: sub div on grid corners convex', '__null__')
      }


      let parts_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '!=', "off_grid");

      for (let part_ of parts_) {

        let subdivs_ = await exec_cluster1_node_q0cbraoo5qe_subdivParts_($p, part_);
        if ($p.terminated) {
          return mfn.getModel();
        }
      }
    }


    async function exec_cluster1_node_ua4tq4tdncr_subdivBlock_($p, block_) {

      let block_id_ = mfn.attrib.Get(block_, 'block_id');

      let block_edges_ = mfn.query.Get('_e', block_);

      let block_verts_ = mfn.query.Get('_v', block_);

      let block_posis_ = mfn.query.Get('ps', block_verts_);

      let cold_edges_ = mfn.query.Filter(mfn.query.Get('_e', block_), ['road', null], '==', "cold");

      let warm_edges_ = mfn.query.Filter(mfn.query.Get('_e', block_), ['road', null], '!=', "cold");

      let road_edges_ = await exec_cluster1_node_ua4tq4tdncr__getEdgesOnParts_($p, block_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      if (ifn.len(road_edges_) == 0) {

        return [];
      }

      let lengths_ = mfn.calc.Length(road_edges_);

      let edge_long_ = ifn.listSort(road_edges_, lengths_)[pythonList(-1, ifn.listSort(road_edges_, lengths_).length)];

      let vec_long_ = mfn.calc.Vector(edge_long_);

      vec_long_ = ifn.vecNorm(vec_long_);

      let vec_long_perp_ = [-vec_long_[pythonList(1, vec_long_.length)], vec_long_[pythonList(0, vec_long_.length)], 0];

      let cut0_ = null;

      let cut1_ = null;

      let posi0_ = null;

      let posi1_ = null;

      let subdivs_ = [];

      for (let i_ of ifn.range(1, ifn.len(road_edges_))) {

        cut0_ = cut1_;

        posi0_ = posi1_;

        let road_edge_ = road_edges_[pythonList(i_, road_edges_.length)];

        let road_idx_ = ifn.listFind(block_edges_, road_edge_);

        let vert_ = mfn.query.Get('_v', road_edge_)[pythonList(0, mfn.query.Get('_v', road_edge_).length)];

        let vert_edges_ = mfn.query.Get('_e', vert_);

        let posi_ = mfn.query.Get('ps', vert_)[pythonList(0, mfn.query.Get('ps', vert_).length)];

        let xyz_ = mfn.attrib.Get(posi_, 'xyz');

        let ang_ = await exec_cluster1_node_ua4tq4tdncr__vertAng_($p, vert_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        let vecs_ = mfn.calc.Vector(vert_edges_);

        vecs_ = ifn.vecNorm(vecs_);

        let vec0_ = ifn.vecRev(vecs_[pythonList(0, vecs_.length)]);

        let vec1_ = vecs_[pythonList(1, vecs_.length)];

        ang_ = ifn.radToDeg(ifn.vecAng2(vec1_, vec0_, [0, 0, 1]));

        if (ang_ > 250) {

          let vec_a_ = ifn.vecRev(ifn.vecSetLen(vec0_, 1000));

          let ray_a_ = ifn.rayMake(xyz_, vec_a_);

          let isect_a_ = await exec_cluster1_node_ua4tq4tdncr__isect_($p, ray_a_, cold_edges_, subdivs_);
          if ($p.terminated) {
            return mfn.getModel();
          }

          let vec_b_ = ifn.vecRev(ifn.vecSetLen(vec1_, 1000));

          let ray_b_ = ifn.rayMake(xyz_, vec_b_);

          let isect_b_ = await exec_cluster1_node_ua4tq4tdncr__isect_($p, ray_b_, cold_edges_, subdivs_);
          if ($p.terminated) {
            return mfn.getModel();
          }

          if (isect_a_[pythonList(0, isect_a_.length)] == null && isect_b_[pythonList(0, isect_b_.length)] == null) {

            continue;
          } else {
            if (isect_a_[pythonList(0, isect_a_.length)] != null && isect_b_[pythonList(0, isect_b_.length)] == null) {

              posi1_ = mfn.make.Position(isect_a_[pythonList(0, isect_a_.length)]);

              let idx_a_ = ifn.listFind(block_edges_, isect_a_[pythonList(1, isect_a_.length)]);

              cut1_ = [road_idx_, idx_a_];
            } else {
              if (isect_a_[pythonList(0, isect_a_.length)] == null && isect_b_[pythonList(0, isect_b_.length)] != null) {

                posi1_ = mfn.make.Position(isect_b_[pythonList(0, isect_b_.length)]);

                let idx_b_ = ifn.listFind(block_edges_, isect_b_[pythonList(1, isect_b_.length)]);

                cut1_ = [road_idx_, idx_b_];
              } else {

                let idx_a_ = ifn.listFind(block_edges_, isect_a_[pythonList(1, isect_a_.length)]);

                let dist_a_ = ifn.distance(isect_a_[pythonList(0, isect_a_.length)], xyz_);

                let idx_b_ = ifn.listFind(block_edges_, isect_b_[pythonList(1, isect_b_.length)]);

                let dist_b_ = ifn.distance(isect_b_[pythonList(0, isect_b_.length)], xyz_);

                if (dist_a_ < dist_b_) {

                  posi1_ = mfn.make.Position(isect_a_[pythonList(0, isect_a_.length)]);

                  cut1_ = [road_idx_, idx_a_];
                } else {

                  posi1_ = mfn.make.Position(isect_b_[pythonList(0, isect_b_.length)]);

                  cut1_ = [road_idx_, idx_b_];
                }
              }
            }
          }
        } else {
          if (ang_ > 190) {

            let vec_avg_ = ifn.vecSum(vecs_);

            let vec_perp_ = ifn.vecSetLen([-vec_avg_[pythonList(1, vec_avg_.length)], vec_avg_[pythonList(0, vec_avg_.length)], 0], 100);

            let ray_perp_ = ifn.rayMake(xyz_, vec_perp_);

            let isect_perp_ = await exec_cluster1_node_ua4tq4tdncr__isect_($p, ray_perp_, cold_edges_, subdivs_);
            if ($p.terminated) {
              return mfn.getModel();
            }

            let idx_perp_ = ifn.listFind(block_edges_, isect_perp_[pythonList(1, isect_perp_.length)]);

            isect_perp_ = isect_perp_[pythonList(0, isect_perp_.length)];

            if (isect_perp_ != null) {

              posi1_ = mfn.make.Position(isect_perp_);

              cut1_ = [road_idx_, idx_perp_];
            } else {

              continue;
            }
          } else {

            continue;
          }
        }

        let subdiv_ = await exec_cluster1_node_ua4tq4tdncr__makeSubdivs_($p, block_posis_, cut0_, cut1_, posi0_, posi1_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (subdiv_ != null) {

          mfn.list.Add(subdivs_, subdiv_, 'to_end');
        }
      }

      posi0_ = posi1_;

      cut0_ = cut1_;

      cut1_ = null;

      let subdiv_ = await exec_cluster1_node_ua4tq4tdncr__makeSubdivs_($p, block_posis_, cut0_, cut1_, posi0_, posi1_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      if (subdiv_ != null) {

        mfn.list.Add(subdivs_, subdiv_, 'to_end');
      }

      await exec_cluster1_node_ua4tq4tdncr__transferEdgeAttribs_($p, block_edges_, mfn.query.Get('_e', subdivs_));
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_node_ua4tq4tdncr__transferEdgeAttribsBtwTouchingParts_($p, subdivs_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_node_ua4tq4tdncr__copyAttribs_($p, block_, subdivs_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      return subdivs_;
    }


    async function exec_cluster1_node_ua4tq4tdncr__transferEdgeAttribs_($p, from_edges_, to_edges_) {

      for (let to_edge_ of to_edges_) {

        let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

        let from_edge_ = await exec_cluster1_node_ua4tq4tdncr__touchingEdge_($p, from_edges_, cen_);
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


    async function exec_cluster1_node_ua4tq4tdncr__transferEdgeAttribsBtwTouchingParts_($p, parts_) {

      let edges_ = mfn.query.Filter(mfn.query.Get('_e', parts_), ['road', null], '==', null);

      for (let to_edge_ of edges_) {

        if (mfn.attrib.Get(to_edge_, 'road') == undefined) {

          let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

          let idx_ = ifn.listFind(edges_, to_edge_);

          let from_edges_ = ifn.listJoin(edges_.slice(0, idx_), edges_.slice(idx_ + 1));

          let from_edge_ = await exec_cluster1_node_ua4tq4tdncr__touchingEdge_($p, from_edges_, cen_);
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


    async function exec_cluster1_node_ua4tq4tdncr__copyAttribs_($p, pgon_from_, pgons_to_) {

      mfn.attrib.Set(pgons_to_, `block_id`, mfn.attrib.Get(pgon_from_, 'block_id'));

      mfn.attrib.Set(pgons_to_, `block_type`, mfn.attrib.Get(pgon_from_, 'block_type'));

      mfn.attrib.Set(pgons_to_, `site`, mfn.attrib.Get(pgon_from_, 'site'));

      mfn.attrib.Set(pgons_to_, `type`, mfn.attrib.Get(pgon_from_, 'type'));

      mfn.attrib.Set(pgons_to_, `class`, mfn.attrib.Get(pgon_from_, 'class'));
    }


    async function exec_cluster1_node_ua4tq4tdncr__touchingEdge_($p, edges_, xyz_) {

      for (let edge_ of edges_) {

        let r_ = mfn.calc.Ray(edge_);

        let d_ = ifn.distance(xyz_, ifn.project(xyz_, r_, 0));

        if (d_ < 1) {

          return edge_;
        }
      }

      return null;
    }


    async function exec_cluster1_node_ua4tq4tdncr__getEdgesOnParts_($p, block_) {

      let edges_ = [[]];

      for (let edge_ of mfn.query.Get('_e', block_)) {

        let edge_road_ = mfn.attrib.Get(edge_, 'road');

        if (edge_road_ == "loc" || edge_road_ == "sec") {

          mfn.list.Add(edges_[pythonList(-1, edges_.length)], edge_, 'to_end');
        } else {

          if (ifn.len(edges_[pythonList(-1, edges_.length)]) != 0) {

            mfn.list.Add(edges_, [], 'to_end');
          }
        }
      }

      if (ifn.len(edges_[pythonList(0, edges_.length)]) == 0) {

        return [];
      }

      if (ifn.len(edges_) == 1) {

        return edges_[pythonList(0, edges_.length)];
      }

      return ifn.listJoin(edges_[pythonList(1, edges_.length)], edges_[pythonList(0, edges_.length)]);
    }


    async function exec_cluster1_node_ua4tq4tdncr__isect_($p, ray_, edges_, subdivs_) {

      let rays_ = mfn.calc.Ray(edges_);

      for (let edge_ of mfn.query.Get('_e', subdivs_)) {

        let edge_ray_ = mfn.calc.Ray(edge_);

        let isect_ = ifn.intersect(ray_, edge_ray_, 0);

        if (isect_ != null) {

          return [null, null];
        }
      }

      let isect_min_ = null;

      let dist_min_ = Infinity;

      let edge_min_ = null;

      for (let edge_ of edges_) {

        let edge_ray_ = mfn.calc.Ray(edge_);

        let isect_ = ifn.intersect(ray_, edge_ray_, 0);

        if (isect_ != null) {

          let dist_ = ifn.distance(isect_, ray_[pythonList(0, ray_.length)]);

          if (dist_ < dist_min_) {

            isect_min_ = isect_;

            dist_min_ = dist_;

            edge_min_ = edge_;
          }
        }
      }

      return [isect_min_, edge_min_];
    }


    async function exec_cluster1_node_ua4tq4tdncr__vertAng_($p, vert_) {

      let edges_ = mfn.query.Get('_e', vert_);

      if (ifn.len(edges_) == 1) {

        return 0;
      }

      let vecs_ = mfn.calc.Vector(edges_);

      let vec0_ = ifn.vecRev(ifn.vecNorm(vecs_[pythonList(0, vecs_.length)]));

      let vec1_ = ifn.vecNorm(vecs_[pythonList(1, vecs_.length)]);

      return ifn.radToDeg(ifn.vecAng2(vec1_, vec0_, [0, 0, 1]));
    }


    async function exec_cluster1_node_ua4tq4tdncr__makeSubdivs_($p, block_posis_, cut0_, cut1_, posi0_, posi1_) {

      let ring_ = [];

      if (cut0_ == null && cut1_ == null) {

        let pgon_ = mfn.make.Polygon(block_posis_);

        mfn.attrib.Set(pgon_, `type`, "off_grid");

        return pgon_;
      } else {
        if (cut0_ == null) {

          ring_ = [posi1_];

          let posis_list_ = await exec_cluster1_node_ua4tq4tdncr__getPosisFromRing_($p, block_posis_, [cut1_[pythonList(1, cut1_.length)] + 1, cut1_[pythonList(0, cut1_.length)] + 1]);
          if ($p.terminated) {
            return mfn.getModel();
          }

          mfn.list.Add(ring_, posis_list_, 'extend_end');
        } else {
          if (cut1_ == null) {

            ring_ = [posi0_];

            let posis_list_ = await exec_cluster1_node_ua4tq4tdncr__getPosisFromRing_($p, block_posis_, [cut0_[pythonList(0, cut0_.length)], cut0_[pythonList(1, cut0_.length)] + 1]);
            if ($p.terminated) {
              return mfn.getModel();
            }

            mfn.list.Add(ring_, posis_list_, 'extend_end');
          } else {

            ring_ = [posi1_];

            let posis_list_ = await exec_cluster1_node_ua4tq4tdncr__getPosisFromRing_($p, block_posis_, [cut1_[pythonList(1, cut1_.length)] + 1, cut0_[pythonList(1, cut0_.length)] + 1]);
            if ($p.terminated) {
              return mfn.getModel();
            }

            mfn.list.Add(ring_, posis_list_, 'extend_end');

            mfn.list.Add(ring_, posi0_, 'extend_end');

            posis_list_ = await exec_cluster1_node_ua4tq4tdncr__getPosisFromRing_($p, block_posis_, [cut0_[pythonList(0, cut0_.length)], cut1_[pythonList(0, cut1_.length)] + 1]);
            if ($p.terminated) {
              return mfn.getModel();
            }

            mfn.list.Add(ring_, posis_list_, 'extend_end');
          }
        }
      }

      if (ifn.len(ring_) < 3) {

        return null;
      }

      let pgon_ = mfn.make.Polygon(ring_);

      mfn.attrib.Set(pgon_, `type`, "off_grid");

      return pgon_;
    }


    async function exec_cluster1_node_ua4tq4tdncr__getPosisFromRing_($p, posis_, idxs_) {

      let idx0_ = idxs_[pythonList(0, idxs_.length)];

      let idx1_ = idxs_[pythonList(1, idxs_.length)];

      if (idx0_ == idx1_) {

        return [];
      }

      let ring_ = [];

      if (idx0_ < idx1_) {

        mfn.list.Add(ring_, posis_.slice(idx0_, idx1_), 'extend_end');
      } else {

        mfn.list.Add(ring_, posis_.slice(idx0_), 'extend_end');

        mfn.list.Add(ring_, posis_.slice(0, idx1_), 'extend_end');
      }

      return ring_;
    }

    async function exec_cluster1_node_ua4tq4tdncr($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: sub div off grid', '__null__')
      }


      let blocks_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', "off_grid");

      for (let block_ of blocks_) {

        let subdivs_ = await exec_cluster1_node_ua4tq4tdncr_subdivBlock_($p, block_);
        if ($p.terminated) {
          return mfn.getModel();
        }
      }

      mfn.edit.Delete([blocks_], 'delete_selected');
    }


    async function exec_cluster1_node_k9k87ylilf_getFrontEdges_($p, cluster_) {

      let pline_ = await exec_cluster1_node_k9k87ylilf__getPerimPlines_($p, cluster_, "loc");
      if ($p.terminated) {
        return mfn.getModel();
      }

      if (ifn.len(pline_) == 0) {

        pline_ = await exec_cluster1_node_k9k87ylilf__getPerimPlines_($p, cluster_, "sec");
        if ($p.terminated) {
          return mfn.getModel();
        }
      }

      if (ifn.len(pline_) == 0) {

        return [[], []];
      }

      let paths_loc_ = [];

      let paths_sec_ = [];

      let pline_edges_ = mfn.query.Get('_e', pline_[pythonList(0, pline_.length)]);

      let posis_loc_ = mfn.query.Get('ps', pline_[pythonList(0, pline_.length)]);

      let lengths_ = mfn.calc.Length(pline_edges_);

      let edge_long_ = ifn.listSort(pline_edges_, lengths_)[pythonList(-1, ifn.listSort(pline_edges_, lengths_).length)];

      let vec_long_ = mfn.calc.Vector(edge_long_);

      vec_long_ = ifn.vecNorm(vec_long_);

      let vec_long_perp_ = [-vec_long_[pythonList(1, vec_long_.length)], vec_long_[pythonList(0, vec_long_.length)], 0];

      let front_edges_ = [edge_long_];

      let idx_ = ifn.listFind(pline_edges_, edge_long_);

      let prev_vec_ = vec_long_;

      for (let i_ of ifn.range(idx_ + 1, ifn.len(pline_edges_), 1)) {

        let edge_ = pline_edges_[pythonList(i_, pline_edges_.length)];

        let vec_ = mfn.calc.Vector(edge_);

        vec_ = ifn.vecNorm(vec_);

        let vec_perp_ = [-vec_[pythonList(1, vec_.length)], vec_[pythonList(0, vec_.length)], 0];

        if (ifn.vecDot(vec_long_, vec_) > 0.8) {

          mfn.list.Add(front_edges_, edge_, 'to_end');

          prev_vec_ = vec_;
        } else {

          break;
        }
      }

      prev_vec_ = vec_long_;

      for (let i_ of ifn.range(idx_ - 1, -1, -1)) {

        let edge_ = pline_edges_[pythonList(i_, pline_edges_.length)];

        let vec_ = mfn.calc.Vector(edge_);

        vec_ = ifn.vecNorm(vec_);

        let vec_perp_ = [-vec_[pythonList(1, vec_.length)], vec_[pythonList(0, vec_.length)], 0];

        if (ifn.vecDot(vec_long_, vec_) > 0.8) {

          mfn.list.Add(front_edges_, edge_, 'to_start');

          prev_vec_ = vec_;
        } else {

          break;
        }
      }

      return front_edges_;
    }


    async function exec_cluster1_node_k9k87ylilf_trimClusterDepth_($p, cluster_, front_edges_) {

      front_edges_ = ifn.listFlat(front_edges_);

      let pline_ = mfn.make.Polyline(mfn.query.Get('ps', front_edges_), 'open');

      let off_ = mfn.poly2d.OffsetMitre(pline_, PART_OG_D_ * 3, 10, 'square_end');

      let remainder_ = mfn.poly2d.Boolean(cluster_, off_, 'difference');

      let remainder_areas_ = mfn.calc.Area(remainder_);

      if (ifn.len(remainder_) == 0 || ifn.sum(remainder_areas_) < (PART_OG_W_ * PART_OG_D_ * 4)) {

        mfn.edit.Delete([pline_, off_, remainder_], 'delete_selected');

        return [];
      }

      let overlap_ = mfn.poly2d.Boolean(cluster_, off_, 'intersect');

      let new_clusters_ = ifn.listJoin(overlap_, remainder_);

      await exec_cluster1_node_k9k87ylilf__copyAttribs_($p, cluster_, new_clusters_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let block_id_ = mfn.attrib.Get(cluster_, 'block_id');

      mfn.attrib.Set(overlap_, `block_id`, block_id_ + "o");

      mfn.attrib.Set(remainder_, `block_id`, block_id_ + "r");

      await exec_cluster1_node_k9k87ylilf__transferEdgeAttribsTouching_($p, mfn.query.Get('_e', cluster_), mfn.query.Get('_e', new_clusters_));
      if ($p.terminated) {
        return mfn.getModel();
      }

      mfn.edit.Delete([pline_, off_, cluster_], 'delete_selected');

      return remainder_;
    }


    async function exec_cluster1_node_k9k87ylilf__getPerimPlines_($p, site_, road_descr_) {

      let posis_ = [];

      for (let edge_ of mfn.query.Get('_e', site_)) {

        let edge_road_ = mfn.attrib.Get(edge_, 'road');

        if (edge_road_ == road_descr_) {

          let start_end_ = mfn.query.Get('ps', edge_);

          if (ifn.len(posis_) == 0 || start_end_[pythonList(0, start_end_.length)] != posis_[pythonList(-1, posis_.length)][pythonList(-1, posis_[pythonList(-1, posis_.length)].length)]) {

            mfn.list.Add(posis_, start_end_, 'to_end');
          } else {

            mfn.list.Add(posis_[pythonList(-1, posis_.length)], start_end_[pythonList(1, start_end_.length)], 'to_end');
          }
        }
      }

      if (ifn.len(posis_) == 0) {

        return [];
      }

      if (ifn.len(posis_) > 1 && posis_[pythonList(-1, posis_.length)][pythonList(-1, posis_[pythonList(-1, posis_.length)].length)] == posis_[pythonList(0, posis_.length)][pythonList(0, posis_[pythonList(0, posis_.length)].length)]) {

        let first_list_ = ifn.listJoin(posis_[pythonList(-1, posis_.length)], posis_[pythonList(0, posis_.length)].slice(1));

        posis_[pythonList(0, posis_.length)] = first_list_;

        posis_ = posis_.slice(0, -1);
      }

      let site_plines_ = mfn.make.Polyline(posis_, 'open');

      return site_plines_;
    }


    async function exec_cluster1_node_k9k87ylilf__transferEdgeAttribsTouching_($p, from_edges_, to_edges_) {

      for (let to_edge_ of to_edges_) {

        let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

        let from_edge_ = await exec_cluster1_node_k9k87ylilf__touchingEdge_($p, from_edges_, cen_);
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


    async function exec_cluster1_node_k9k87ylilf__copyAttribs_($p, pgon_from_, pgons_to_) {

      let block_id_ = mfn.attrib.Get(pgon_from_, 'block_id');

      let block_type_ = mfn.attrib.Get(pgon_from_, 'block_type');

      let site_ = mfn.attrib.Get(pgon_from_, 'site');

      mfn.attrib.Set(pgons_to_, `block_id`, mfn.attrib.Get(pgon_from_, 'block_id'));

      mfn.attrib.Set(pgons_to_, `block_type`, mfn.attrib.Get(pgon_from_, 'block_type'));

      mfn.attrib.Set(pgons_to_, `site`, mfn.attrib.Get(pgon_from_, 'site'));

      mfn.attrib.Set(pgons_to_, `class`, mfn.attrib.Get(pgon_from_, 'class'));

      mfn.attrib.Set(pgons_to_, `type`, mfn.attrib.Get(pgon_from_, 'type'));
    }


    async function exec_cluster1_node_k9k87ylilf__touchingEdge_($p, edges_, xyz_) {

      for (let edge_ of edges_) {

        let r_ = mfn.calc.Ray(edge_);

        let d_ = ifn.distance(xyz_, ifn.project(xyz_, r_, 0));

        if (d_ < 1) {

          return edge_;
        }
      }

      return null;
    }

    async function exec_cluster1_node_k9k87ylilf($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: split off grid too deep', '__null__')
      }


      let clusters_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', "off_grid");

      clusters_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', "off_grid");

      let new_clusters1_ = [];

      for (let cluster_ of clusters_) {

        let front_edges_ = await exec_cluster1_node_k9k87ylilf_getFrontEdges_($p, cluster_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        let remainder_ = await exec_cluster1_node_k9k87ylilf_trimClusterDepth_($p, cluster_, front_edges_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.list.Add(new_clusters1_, remainder_, 'extend_end');
      }

      let new_clusters2_ = [];

      for (let cluster_ of new_clusters1_) {

        let front_edges_ = await exec_cluster1_node_k9k87ylilf_getFrontEdges_($p, cluster_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        let remainder_ = await exec_cluster1_node_k9k87ylilf_trimClusterDepth_($p, cluster_, front_edges_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.list.Add(new_clusters2_, remainder_, 'extend_end');
      }

      for (let cluster_ of new_clusters2_) {

        let front_edges_ = await exec_cluster1_node_k9k87ylilf_getFrontEdges_($p, cluster_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        let remainder_ = await exec_cluster1_node_k9k87ylilf_trimClusterDepth_($p, cluster_, front_edges_);
        if ($p.terminated) {
          return mfn.getModel();
        }
      }

      let clusters3_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', "off_grid");

      for (let cluster_ of clusters3_) {

        let edges_ = mfn.query.Get('_e', cluster_);

        if (ifn.len(mfn.query.Filter(edges_, ['road', null], '==', "loc")) > 0) {

          continue;
        }

        if (ifn.len(mfn.query.Filter(edges_, ['road', null], '==', "sec")) > 0) {

          continue;
        }

        mfn.attrib.Set(cluster_, `type`, "too_deep");

        mfn.attrib.Set(cluster_, `class`, "leftover");
      }
    }


    async function exec_cluster1_node_sqstwmuvkl8_processCluster_($p, cluster_) {

      let new_parts_list_ = [];

      let area_ = mfn.calc.Area(cluster_);

      let front_edges_ = await exec_cluster1_node_sqstwmuvkl8__getFrontEdges_($p, cluster_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      if (ifn.len(front_edges_) == 0) {

        return [];
      }

      if (area_ > (PART_OG_D_ * PART_OG_W_ * 2)) {

        let new_parts_lists_ = await exec_cluster1_node_sqstwmuvkl8__subdivFrontEdges_($p, cluster_, front_edges_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        new_parts_list_ = ifn.listFlat(new_parts_lists_);

        mfn.edit.Delete(cluster_, 'delete_selected');
      } else {

        let cluster_edges_ = mfn.query.Get('_e', cluster_);

        let front_edge_ = await exec_cluster1_node_sqstwmuvkl8__getLongestEdge_($p, front_edges_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        let cen_ = mfn.calc.Centroid(front_edge_, 'ps_average');

        let edge_ = await exec_cluster1_node_sqstwmuvkl8__touchingEdge_($p, cluster_edges_, cen_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        let front_idx_ = ifn.listFind(cluster_edges_, edge_);

        if (front_idx_ != null) {

          mfn.edit.Shift(cluster_, front_idx_);
        }

        mfn.attrib.Set(cluster_, `type`, "off_grid1");

        mfn.attrib.Set(cluster_, `entr`, cen_);

        new_parts_list_ = [cluster_];
      }

      mfn.attrib.Set(new_parts_list_, `class`, "part");

      return new_parts_list_;
    }


    async function exec_cluster1_node_sqstwmuvkl8__getLongestEdge_($p, edges_) {

      if (ifn.len(edges_) == 1) {

        return edges_[pythonList(0, edges_.length)];
      }

      let lengths_ = mfn.calc.Length(edges_);

      let long_edge_ = ifn.listSort(edges_, lengths_)[pythonList(-1, ifn.listSort(edges_, lengths_).length)];

      return long_edge_;
    }


    async function exec_cluster1_node_sqstwmuvkl8__getFrontEdges_($p, cluster_) {

      let pline_ = await exec_cluster1_node_sqstwmuvkl8__getPerimPlines_($p, cluster_, "loc");
      if ($p.terminated) {
        return mfn.getModel();
      }

      if (ifn.len(pline_) == 0) {

        pline_ = await exec_cluster1_node_sqstwmuvkl8__getPerimPlines_($p, cluster_, "sec");
        if ($p.terminated) {
          return mfn.getModel();
        }
      }

      if (ifn.len(pline_) == 0) {

        pline_ = await exec_cluster1_node_sqstwmuvkl8__getPerimPlines_($p, cluster_, "art");
        if ($p.terminated) {
          return mfn.getModel();
        }
      }

      if (ifn.len(pline_) == 0) {

        return [];
      }

      let paths_loc_ = [];

      let paths_sec_ = [];

      let pline_edges_ = mfn.query.Get('_e', pline_[pythonList(0, pline_.length)]);

      let posis_loc_ = mfn.query.Get('ps', pline_[pythonList(0, pline_.length)]);

      let lengths_ = mfn.calc.Length(pline_edges_);

      let edge_long_ = ifn.listSort(pline_edges_, lengths_)[pythonList(-1, ifn.listSort(pline_edges_, lengths_).length)];

      let vec_long_ = mfn.calc.Vector(edge_long_);

      vec_long_ = ifn.vecNorm(vec_long_);

      let vec_long_perp_ = [-vec_long_[pythonList(1, vec_long_.length)], vec_long_[pythonList(0, vec_long_.length)], 0];

      let front_edges_ = [edge_long_];

      let idx_ = ifn.listFind(pline_edges_, edge_long_);

      let prev_vec_ = vec_long_;

      for (let i_ of ifn.range(idx_ + 1, ifn.len(pline_edges_), 1)) {

        let edge_ = pline_edges_[pythonList(i_, pline_edges_.length)];

        let vec_ = mfn.calc.Vector(edge_);

        vec_ = ifn.vecNorm(vec_);

        let vec_perp_ = [-vec_[pythonList(1, vec_.length)], vec_[pythonList(0, vec_.length)], 0];

        if (ifn.vecDot(vec_long_, vec_) > 0.8) {

          mfn.list.Add(front_edges_, edge_, 'to_end');

          prev_vec_ = vec_;
        } else {

          break;
        }
      }

      prev_vec_ = vec_long_;

      for (let i_ of ifn.range(idx_ - 1, -1, -1)) {

        let edge_ = pline_edges_[pythonList(i_, pline_edges_.length)];

        let vec_ = mfn.calc.Vector(edge_);

        vec_ = ifn.vecNorm(vec_);

        let vec_perp_ = [-vec_[pythonList(1, vec_.length)], vec_[pythonList(0, vec_.length)], 0];

        if (ifn.vecDot(vec_long_, vec_) > 0.8) {

          mfn.list.Add(front_edges_, edge_, 'to_start');

          prev_vec_ = vec_;
        } else {

          break;
        }
      }

      return front_edges_;
    }


    async function exec_cluster1_node_sqstwmuvkl8__subdivFrontEdges_($p, cluster_, front_edges_) {

      let cluster_edges_ = mfn.query.Get('_e', cluster_);

      let lengths_ = mfn.calc.Length(front_edges_);

      let width_tot_ = ifn.sum(lengths_);

      let num_og_ = ifn.round(width_tot_ / PART_OG_W_);

      let width_og_ = width_tot_ / num_og_;

      let half_width_og_ = width_og_ / 2;

      let cut_rays_ = [];

      let entrs_ = [];

      let start_offset_ = 0;

      for (let edge_ of front_edges_) {

        let xyzs_ = mfn.attrib.Get(mfn.query.Get('ps', edge_), 'xyz');

        let vec_ = mfn.calc.Vector(edge_);

        let vec_perp_ = ifn.vecSetLen([-vec_[pythonList(1, vec_.length)], vec_[pythonList(0, vec_.length)], 0], 500);

        let edge_length_ = ifn.vecLen(vec_);

        let div_length_ = edge_length_ - start_offset_;

        let num_div_ = ifn.floor((div_length_ + 0.001) / half_width_og_) + 1;

        let ori_ = ifn.vecAdd(xyzs_[pythonList(0, xyzs_.length)], ifn.vecSetLen(vec_, start_offset_));

        for (let i_ of ifn.range(num_div_)) {

          let cut_ray_ = ifn.rayMake(ori_, vec_perp_);

          if (edge_ == front_edges_[pythonList(0, front_edges_.length)] && i_ == 0) {

            cut_ray_[pythonList(0, cut_ray_.length)] = ifn.vecAdd(cut_ray_[pythonList(0, cut_ray_.length)], ifn.vecSetLen(vec_, -100));
          } else {
            if (edge_ == front_edges_[pythonList(-1, front_edges_.length)] && i_ == num_div_ - 1) {

              cut_ray_[pythonList(0, cut_ray_.length)] = ifn.vecAdd(cut_ray_[pythonList(0, cut_ray_.length)], ifn.vecSetLen(vec_, 100));
            }
          }

          mfn.list.Add(cut_rays_, cut_ray_, 'to_end');

          ori_ = ifn.vecAdd(ori_, ifn.vecSetLen(vec_, half_width_og_));
        }

        let div_sum_ = (num_div_ - 1) * half_width_og_;

        start_offset_ = half_width_og_ - (div_length_ - div_sum_);
      }

      let new_parts_list_ = [];

      let white_parts_ = [];

      for (let i_ of ifn.range(0, ifn.len(cut_rays_) - 2, 2)) {

        let ray0_ = ifn.rayLMove(cut_rays_[pythonList(i_, cut_rays_.length)], -10);

        let ray1_ = ifn.rayLMove(cut_rays_[pythonList(i_ + 2, cut_rays_.length)], -10);

        let entr_ = cut_rays_[pythonList(i_ + 1, cut_rays_.length)][pythonList(0, cut_rays_[pythonList(i_ + 1, cut_rays_.length)].length)];

        let pgon1_ = await exec_cluster1_node_sqstwmuvkl8__makePgonFromRays_($p, ray0_, ray1_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        let og_ = mfn.poly2d.Boolean(pgon1_, cluster_, 'intersect');

        let area_og_ = mfn.calc.Area(og_);

        let two_parts_ = null;

        if (area_og_ > (PART_OG_D_ * PART_OG_W_ * 1.5)) {

          let vec_ = ifn.vecSetLen(ifn.vecFromTo(ray0_[pythonList(0, ray0_.length)], ray1_[pythonList(0, ray1_.length)]), 1);

          ray0_[pythonList(0, ray0_.length)] = ifn.vecSub(ray0_[pythonList(0, ray0_.length)], vec_);

          ray1_[pythonList(0, ray1_.length)] = ifn.vecAdd(ray1_[pythonList(0, ray1_.length)], vec_);

          ray0_[pythonList(1, ray0_.length)] = ifn.vecSetLen(ray0_[pythonList(1, ray0_.length)], PART_OG_D_ + 10);

          ray1_[pythonList(1, ray1_.length)] = ifn.vecSetLen(ray1_[pythonList(1, ray1_.length)], PART_OG_D_ + 10);

          let pgon2_ = await exec_cluster1_node_sqstwmuvkl8__makePgonFromRays_($p, ray0_, ray1_);
          if ($p.terminated) {
            return mfn.getModel();
          }

          let og1_ = mfn.poly2d.Boolean(og_, pgon2_, 'intersect');

          let og2_ = mfn.poly2d.Boolean(og_, pgon2_, 'difference');

          mfn.edit.Delete([pgon1_, pgon2_, og_], 'delete_selected');

          mfn.attrib.Set(og1_, `entr`, entr_);

          if (ifn.len(og1_) > 0 && ifn.len(og2_) > 0) {

            let cluster_id_ = mfn.attrib.Get(cluster_, 'block_id') + '_' + i_;

            mfn.attrib.Set(og1_, `cluster_id`, cluster_id_);

            mfn.attrib.Set(og2_, `cluster_id`, cluster_id_);
          }

          two_parts_ = [og1_, og2_];
        } else {

          mfn.attrib.Set(og_, `entr`, entr_);

          two_parts_ = [og_, []];

          mfn.edit.Delete(pgon1_, 'delete_selected');
        }

        let whites_ = await exec_cluster1_node_sqstwmuvkl8__filterAndSortParts_($p, two_parts_, entr_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.list.Add(white_parts_, whites_, 'extend_end');

        let part0_ = two_parts_[pythonList(0, two_parts_.length)][pythonList(0, two_parts_[pythonList(0, two_parts_.length)].length)];

        let part0_edges_ = mfn.query.Get('_e', part0_);

        let edge_ = await exec_cluster1_node_sqstwmuvkl8__touchingEdge_($p, part0_edges_, mfn.attrib.Get(part0_, 'entr'));
        if ($p.terminated) {
          return mfn.getModel();
        }

        let front_idx_ = ifn.listFind(part0_edges_, edge_);

        if (front_idx_ != null) {

          mfn.edit.Shift(part0_, front_idx_);
        }

        mfn.list.Add(new_parts_list_, two_parts_, 'to_end');
      }

      await exec_cluster1_node_sqstwmuvkl8__attribs_($p, cluster_, new_parts_list_, white_parts_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      return new_parts_list_;
    }


    async function exec_cluster1_node_sqstwmuvkl8__attribs_($p, cluster_, new_parts_list_, white_parts_) {

      for (let two_parts_ of new_parts_list_) {

        await exec_cluster1_node_sqstwmuvkl8__copyAttribs_($p, cluster_, ifn.listFlat(two_parts_));
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.attrib.Set(two_parts_[pythonList(0, two_parts_.length)], `type`, "off_grid1");

        mfn.attrib.Set(two_parts_[pythonList(1, two_parts_.length)], `type`, "off_grid2");
      }

      await exec_cluster1_node_sqstwmuvkl8__copyAttribs_($p, cluster_, white_parts_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      mfn.attrib.Set(white_parts_, `type`, "white");

      let all_parts_ = ifn.listFlat([new_parts_list_, white_parts_]);

      mfn.attrib.Set(all_parts_, `class`, "part");

      await exec_cluster1_node_sqstwmuvkl8__transferEdgeAttribsTouching_($p, mfn.query.Get('_e', cluster_), mfn.query.Get('_e', all_parts_));
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_node_sqstwmuvkl8__transferEdgeAttribsBtwTouchingParts_($p, all_parts_);
      if ($p.terminated) {
        return mfn.getModel();
      }
    }


    async function exec_cluster1_node_sqstwmuvkl8__makePgonFromRays_($p, ray0_, ray1_) {

      let posis_ = mfn.make.Position([ray0_[pythonList(0, ray0_.length)], ray1_[pythonList(0, ray1_.length)], ifn.vecAdd(ray1_[pythonList(0, ray1_.length)], ray1_[pythonList(1, ray1_.length)]), ifn.vecAdd(ray0_[pythonList(0, ray0_.length)], ray0_[pythonList(1, ray0_.length)])]);

      let pgon_ = mfn.make.Polygon(posis_);

      return pgon_;
    }


    async function exec_cluster1_node_sqstwmuvkl8__filterAndSortParts_($p, two_parts_, entr_) {

      let whites_ = [];

      for (let i_ of [0, 1]) {

        let parts_ = two_parts_[pythonList(i_, two_parts_.length)];

        if (ifn.len(parts_) < 2) {

          continue;
        }

        let part_dists_ = [];

        let part_sort_ = [];

        for (let part_ of parts_) {

          let cen_ = mfn.calc.Centroid(part_, 'ps_average');

          let dist_ = ifn.distance(entr_, cen_);

          mfn.list.Add(part_dists_, [part_, dist_], 'to_end');

          mfn.list.Add(part_sort_, dist_, 'to_end');
        }

        let sorted_ = ifn.listSort(part_dists_, part_sort_);

        for (let part_dist_ of sorted_.slice(1)) {

          mfn.list.Add(whites_, part_dist_[pythonList(0, part_dist_.length)], 'to_end');
        }

        two_parts_[pythonList(i_, two_parts_.length)] = [sorted_[pythonList(0, sorted_.length)][pythonList(0, sorted_[pythonList(0, sorted_.length)].length)]];
      }

      return whites_;
    }


    async function exec_cluster1_node_sqstwmuvkl8__getPerimPlines_($p, site_, road_descr_) {

      let posis_ = [];

      for (let edge_ of mfn.query.Get('_e', site_)) {

        let edge_road_ = mfn.attrib.Get(edge_, 'road');

        if (edge_road_ == road_descr_) {

          let start_end_ = mfn.query.Get('ps', edge_);

          if (ifn.len(posis_) == 0 || start_end_[pythonList(0, start_end_.length)] != posis_[pythonList(-1, posis_.length)][pythonList(-1, posis_[pythonList(-1, posis_.length)].length)]) {

            mfn.list.Add(posis_, start_end_, 'to_end');
          } else {

            mfn.list.Add(posis_[pythonList(-1, posis_.length)], start_end_[pythonList(1, start_end_.length)], 'to_end');
          }
        }
      }

      if (ifn.len(posis_) == 0) {

        return [];
      }

      if (ifn.len(posis_) > 1 && posis_[pythonList(-1, posis_.length)][pythonList(-1, posis_[pythonList(-1, posis_.length)].length)] == posis_[pythonList(0, posis_.length)][pythonList(0, posis_[pythonList(0, posis_.length)].length)]) {

        let first_list_ = ifn.listJoin(posis_[pythonList(-1, posis_.length)], posis_[pythonList(0, posis_.length)].slice(1));

        posis_[pythonList(0, posis_.length)] = first_list_;

        posis_ = posis_.slice(0, -1);
      }

      let site_plines_ = mfn.make.Polyline(posis_, 'open');

      return site_plines_;
    }


    async function exec_cluster1_node_sqstwmuvkl8__transferEdgeAttribsTouching_($p, from_edges_, to_edges_) {

      for (let to_edge_ of to_edges_) {

        let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

        let from_edge_ = await exec_cluster1_node_sqstwmuvkl8__touchingEdge_($p, from_edges_, cen_);
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


    async function exec_cluster1_node_sqstwmuvkl8__copyAttribs_($p, pgon_from_, pgons_to_) {

      let block_id_ = mfn.attrib.Get(pgon_from_, 'block_id');

      let block_type_ = mfn.attrib.Get(pgon_from_, 'block_type');

      let site_ = mfn.attrib.Get(pgon_from_, 'site');

      for (let pgon_to_ of pgons_to_) {

        mfn.attrib.Set(pgon_to_, `block_id`, block_id_);

        mfn.attrib.Set(pgon_to_, `block_type`, block_type_);

        mfn.attrib.Set(pgon_to_, `site`, site_);
      }
    }


    async function exec_cluster1_node_sqstwmuvkl8__touchingEdge_($p, edges_, xyz_) {

      for (let edge_ of edges_) {

        let r_ = mfn.calc.Ray(edge_);

        let d_ = ifn.distance(xyz_, ifn.project(xyz_, r_, 0));

        if (d_ < 1) {

          return edge_;
        }
      }

      return null;
    }


    async function exec_cluster1_node_sqstwmuvkl8__project_($p, ray_, rays_, max_dist_) {

      let min_dist_ = max_dist_;

      let min_i_ = null;

      let min_isect_ = null;

      for (let i_ of ifn.range(ifn.len(rays_))) {

        let isect_ = ifn.project(ray_[pythonList(0, ray_.length)], rays_[pythonList(i_, rays_.length)], 0);

        if (isect_ != null) {

          if (ifn.vecDot(ifn.vecNorm(ifn.vecFromTo(ray_[pythonList(0, ray_.length)], isect_)), ray_[pythonList(1, ray_.length)]) > 0.7) {

            let dist_ = ifn.distance(ray_[pythonList(0, ray_.length)], isect_);

            if (dist_ < min_dist_) {

              min_dist_ = dist_;

              min_i_ = i_;

              min_isect_ = isect_;
            }
          }
        }
      }

      if (min_i_ == null) {

        return null;
      }

      return [min_i_, min_isect_];
    }


    async function exec_cluster1_node_sqstwmuvkl8__isect_($p, ray_, rays_) {

      for (let i_ of ifn.range(ifn.len(rays_))) {

        let isect_ = ifn.intersect(ray_, rays_[pythonList(i_, rays_.length)], 0);

        if (isect_ != null) {

          return [i_, isect_];
        }
      }

      return null;
    }


    async function exec_cluster1_node_sqstwmuvkl8__vertAng_($p, vert_) {

      let edges_ = mfn.query.Get('_e', vert_);

      if (ifn.len(edges_) == 1) {

        return 0;
      }

      let vecs_ = mfn.calc.Vector(edges_);

      let vec0_ = ifn.vecRev(ifn.vecNorm(vecs_[pythonList(0, vecs_.length)]));

      let vec1_ = ifn.vecNorm(vecs_[pythonList(1, vecs_.length)]);

      return ifn.radToDeg(ifn.vecAng2(vec1_, vec0_, [0, 0, 1]));
    }


    async function exec_cluster1_node_sqstwmuvkl8__transferEdgeAttribsBtwTouchingParts_($p, parts_) {

      let edges_ = mfn.query.Filter(mfn.query.Get('_e', parts_), ['road', null], '==', null);

      for (let to_edge_ of edges_) {

        if (mfn.attrib.Get(to_edge_, 'road') == undefined) {

          let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

          let idx_ = ifn.listFind(edges_, to_edge_);

          let from_edges_ = ifn.listJoin(edges_.slice(0, idx_), edges_.slice(idx_ + 1));

          let from_edge_ = await exec_cluster1_node_sqstwmuvkl8__touchingEdge_($p, from_edges_, cen_);
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

    async function exec_cluster1_node_sqstwmuvkl8($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: create off grid cold', '__null__')
      }


      let clusters_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', "off_grid");

      let clusterps_ = mfn.query.Get('ps', clusters_);

      let newclusters_ = mfn.make.Clone(clusters_);

      for (let cluster_ of newclusters_) {
        printFunc($p.console, 'Executing For-each: cluster = ' + (cluster_), '__null__');

        let new_parts_list_ = await exec_cluster1_node_sqstwmuvkl8_processCluster_($p, cluster_);
        if ($p.terminated) {
          return mfn.getModel();
        }
      }
    }


    async function exec_cluster1_node_cl6rfx4972q_mergeParts_($p, parts_, expected_areas_) {

      let new_parts_ = [];

      for (let part_ of parts_) {

        let exists_ = mfn.query.Type(part_, 'exists');

        if (!exists_) {

          continue;
        }

        let area_ = mfn.calc.Area(part_);

        if (area_ < 1) {

          continue;
        }

        let type_ = mfn.attrib.Get(part_, 'type');

        let expected_area_ = expected_areas_[pythonList(type_, expected_areas_.length)];

        area_ = mfn.calc.Area(part_);

        if (expected_area_ != undefined && area_ < (expected_area_ * 0.6)) {

          let new_part_ = await exec_cluster1_node_cl6rfx4972q__mergeWithNeighbour_($p, part_, parts_);
          if ($p.terminated) {
            return mfn.getModel();
          }

          if (new_part_ != null) {

            mfn.visualize.Color(new_part_, [1, 1, 0]);

            mfn.list.Add(new_parts_, new_part_, 'to_end');
          }
        }
      }

      return new_parts_;
    }


    async function exec_cluster1_node_cl6rfx4972q_getLongestEdge_($p, edges_) {

      if (ifn.len(edges_) == 1) {

        return edges_[pythonList(0, edges_.length)];
      }

      let lengths_ = mfn.calc.Length(edges_);

      let long_edge_ = ifn.listSort(edges_, lengths_)[pythonList(-1, ifn.listSort(edges_, lengths_).length)];

      return long_edge_;
    }


    async function exec_cluster1_node_cl6rfx4972q__mergeWithNeighbour_($p, part_, parts_) {

      let exists_ = mfn.query.Type(part_, 'exists');

      if (!exists_) {

        return null;
      }

      let pairs_ = [];

      let pairs_sorter_ = [];

      for (let edge_ of mfn.query.Get('_e', part_)) {

        let nedge_ = await exec_cluster1_node_cl6rfx4972q__getNeighbourEdge_($p, edge_, parts_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (nedge_ != null) {

          let length_ = mfn.calc.Length(edge_);

          mfn.list.Add(pairs_, [edge_, nedge_, length_], 'to_end');

          mfn.list.Add(pairs_sorter_, length_, 'to_end');
        }
      }

      if (ifn.len(pairs_) > 0) {

        let sorted_ = ifn.listSort(pairs_, pairs_sorter_);

        let longest_edge_ = sorted_[pythonList(-1, sorted_.length)];

        let edge_posis_ = mfn.query.Get('ps', longest_edge_[pythonList(0, longest_edge_.length)]);

        let npart_ = mfn.query.Get('pg', longest_edge_[pythonList(1, longest_edge_.length)])[pythonList(0, mfn.query.Get('pg', longest_edge_[pythonList(1, longest_edge_.length)]).length)];

        let posi_a_ = edge_posis_[pythonList(0, edge_posis_.length)];

        let posi_b_ = edge_posis_[pythonList(1, edge_posis_.length)];

        let posis0_ = mfn.query.Get('ps', part_);

        let len_posis0_ = ifn.len(posis0_);

        posis0_ = ifn.listJoin(posis0_, posis0_);

        let posis1_ = mfn.query.Get('ps', npart_);

        let len_posis1_ = ifn.len(posis1_);

        posis1_ = ifn.listJoin(posis1_, posis1_);

        let idx0_start_ = ifn.listFind(posis0_, posi_b_);

        let idx0_end_ = ifn.listFind(posis0_, posi_a_);

        if (idx0_end_ < idx0_start_) {

          idx0_end_ = idx0_end_ + len_posis0_;
        }

        let idx1_start_ = ifn.listFind(posis1_, posi_a_);

        let idx1_end_ = ifn.listFind(posis1_, posi_b_);

        if (idx1_end_ < idx1_start_) {

          idx1_end_ = idx1_end_ + len_posis1_;
        }

        let posis0_c_ = posis0_.slice(idx0_start_, idx0_end_);

        let posis1_c_ = posis1_.slice(idx1_start_, idx1_end_);

        let chk0_ = [idx0_start_, idx0_end_];

        let chk1_ = [idx1_start_, idx1_end_];

        let posis_ = ifn.listJoin(posis0_.slice(idx0_start_, idx0_end_), posis1_.slice(idx1_start_, idx1_end_));

        await exec_cluster1_node_cl6rfx4972q__removeCracks_($p, posis_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        let new_part_ = mfn.make.Polygon(posis_);

        await exec_cluster1_node_cl6rfx4972q__copyAttribs_($p, npart_, new_part_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        let old_parts_ = ifn.listFlat([part_, npart_]);

        await exec_cluster1_node_cl6rfx4972q__transferEdgeAttribsTouching_($p, mfn.query.Get('_e', old_parts_), mfn.query.Get('_e', new_part_));
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.edit.Delete(old_parts_, 'delete_selected');

        return new_part_;
      }
    }


    async function exec_cluster1_node_cl6rfx4972q__removeCracks_($p, posis_) {

      let cull_ = null;

      let len_posis_ = ifn.len(posis_);

      for (let i_ of ifn.range(ifn.len(posis_))) {

        if (posis_[pythonList(i_, posis_.length)] == posis_[pythonList((i_ + 2) % len_posis_, posis_.length)]) {

          cull_ = [(i_ + 1) % len_posis_, (i_ + 2) % len_posis_];

          break;
        }
      }

      if (cull_ != null) {

        mfn.list.Remove(posis_, cull_[pythonList(1, cull_.length)], 'index');

        mfn.list.Remove(posis_, cull_[pythonList(0, cull_.length)], 'index');

        await exec_cluster1_node_cl6rfx4972q__removeCracks_($p, posis_);
        if ($p.terminated) {
          return mfn.getModel();
        }
      }
    }


    async function exec_cluster1_node_cl6rfx4972q__getNeighbourEdge_($p, edge_, parts_) {

      let posis_ = mfn.query.Get('ps', edge_);

      let edges_ = mfn.query.Get('_e', posis_);

      for (let nedge_ of edges_) {

        if (nedge_ == edge_) {

          continue;
        }

        if (!ifn.listHas(parts_, mfn.query.Get('pg', nedge_)[pythonList(0, mfn.query.Get('pg', nedge_).length)])) {

          continue;
        }

        let nposis_ = mfn.query.Get('ps', nedge_);

        if (nposis_[pythonList(0, nposis_.length)] == posis_[pythonList(1, posis_.length)] && nposis_[pythonList(1, nposis_.length)] == posis_[pythonList(0, posis_.length)]) {

          return nedge_;
        }
      }

      return null;
    }


    async function exec_cluster1_node_cl6rfx4972q__copyAttribs_($p, pgon_from_, pgon_to_) {

      mfn.attrib.Set(pgon_to_, `block_id`, mfn.attrib.Get(pgon_from_, 'block_id'));

      mfn.attrib.Set(pgon_to_, `block_type`, mfn.attrib.Get(pgon_from_, 'block_type'));

      mfn.attrib.Set(pgon_to_, `site`, mfn.attrib.Get(pgon_from_, 'site'));

      mfn.attrib.Set(pgon_to_, `type`, mfn.attrib.Get(pgon_from_, 'type'));

      mfn.attrib.Set(pgon_to_, `class`, mfn.attrib.Get(pgon_from_, 'class'));
    }


    async function exec_cluster1_node_cl6rfx4972q__transferEdgeAttribsTouching_($p, from_edges_, to_edges_) {

      for (let to_edge_ of to_edges_) {

        let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

        let from_edge_ = await exec_cluster1_node_cl6rfx4972q__touchingEdge_($p, from_edges_, cen_);
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


    async function exec_cluster1_node_cl6rfx4972q__touchingEdge_($p, edges_, xyz_) {

      for (let edge_ of edges_) {

        let r_ = mfn.calc.Ray(edge_);

        let d_ = ifn.distance(xyz_, ifn.project(xyz_, r_, 0));

        if (d_ < 1) {

          return edge_;
        }
      }

      return null;
    }


    async function exec_cluster1_node_cl6rfx4972q__cleanPgonsAng_($p, pgons_) {

      for (let pgon_ of pgons_) {

        let var_ = mfn.edit.Weld(pgon_, 'break_weld');

        let del_posis_ = [];

        for (let vert_ of mfn.query.Get('_v', pgon_)) {

          let dot_ = await exec_cluster1_node_cl6rfx4972q__angDot_($p, vert_);
          if ($p.terminated) {
            return mfn.getModel();
          }

          if (ifn.abs(dot_) > 0.9999) {

            mfn.list.Add(del_posis_, mfn.query.Get('ps', vert_), 'to_end');
          }
        }
      }

      return mfn.query.Get('pg', pgons_);
    }


    async function exec_cluster1_node_cl6rfx4972q__angDot_($p, vert_) {

      let edges_ = mfn.query.Get('_e', vert_);

      let vecs_ = mfn.calc.Vector(edges_);

      let vec0_ = ifn.vecNorm(vecs_[pythonList(0, vecs_.length)]);

      let vec1_ = ifn.vecNorm(vecs_[pythonList(1, vecs_.length)]);

      return ifn.vecDot(vec0_, vec1_);
    }

    async function exec_cluster1_node_cl6rfx4972q($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: merge small pgons off_grid1', '__null__')
      }


      let allpgons_ = mfn.query.Get('pg', null);

      let pg_ = mfn.make.Clone(allpgons_);

      mfn.edit.Delete(mfn.query.Get('pl', null), 'delete_selected');

      let expected_areas_ = [];

      expected_areas_["off_grid1"] = PART_OG_D_ * PART_OG_W_;

      let var_ = mfn.edit.Fuse(pg_, 0.01);

      let parts_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', "off_grid1");

      let new_parts1_ = await exec_cluster1_node_cl6rfx4972q_mergeParts_($p, parts_, expected_areas_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      parts_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', "off_grid1");

      let new_parts2_ = await exec_cluster1_node_cl6rfx4972q_mergeParts_($p, parts_, expected_areas_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let new_parts_ = ifn.listJoin(new_parts1_, new_parts2_);

      new_parts_ = mfn.query.Get('pg', new_parts_);

      new_parts_ = await exec_cluster1_node_cl6rfx4972q__cleanPgonsAng_($p, new_parts_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      for (let new_part_ of new_parts_) {

        let edges_ = mfn.query.Get('_e', new_part_);

        let front_edges_ = mfn.query.Filter(edges_, ['road', null], '==', "loc");

        if (ifn.len(front_edges_) == 0) {

          front_edges_ = mfn.query.Filter(edges_, ['road', null], '==', "sec");
        }

        if (ifn.len(front_edges_) == 0) {

          front_edges_ = mfn.query.Filter(edges_, ['road', null], '==', "art");
        }

        if (ifn.len(front_edges_) == 0) {

          continue;
        }

        let front_edge_ = await exec_cluster1_node_cl6rfx4972q_getLongestEdge_($p, front_edges_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.edit.Shift(new_part_, ifn.listFind(edges_, front_edge_));
      }
    }


    async function exec_cluster1_node_zexd2z8qt7($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: get concave corners', '__null__')
      }


      let blocks_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', "block_corner");

      mfn.edit.Delete(blocks_, 'keep_selected');
    }


    async function exec_cluster1_node_agq2g6vtb3_processBlock_($p, block_) {

      let roads_art_ = await exec_cluster1_node_agq2g6vtb3__getSitePlines_($p, block_, "road_art");
      if ($p.terminated) {
        return mfn.getModel();
      }

      let roads_sec_ = await exec_cluster1_node_agq2g6vtb3__getSitePlines_($p, block_, "road_sec");
      if ($p.terminated) {
        return mfn.getModel();
      }

      let result_ = await exec_cluster1_node_agq2g6vtb3__createOnGridStrips_($p, block_, roads_art_, roads_sec_, []);
      if ($p.terminated) {
        return mfn.getModel();
      }

      return result_;
    }


    async function exec_cluster1_node_agq2g6vtb3__createOnGridStrips_($p, block_, roads_art_, roads_sec_, roads_loc_) {

      let off_art_ = mfn.poly2d.OffsetMitre(roads_art_, PART_ART_D_, 100, 'square_end');

      let off_sec_ = mfn.poly2d.OffsetMitre(roads_sec_, PART_SEC_D_, 100, 'square_end');

      let off_loc_ = mfn.poly2d.OffsetMitre(roads_loc_, PART_LOC_D_, 100, 'square_end');

      let off_grids_ = mfn.poly2d.Boolean(block_, [off_art_, off_sec_, off_loc_], 'difference');

      let on_arts1_ = mfn.poly2d.Boolean(off_art_, block_, 'intersect');

      let on_secs1_ = mfn.poly2d.Boolean(off_sec_, block_, 'intersect');

      let on_locs1_ = mfn.poly2d.Boolean(off_loc_, block_, 'intersect');

      let on_arts2_ = mfn.poly2d.Boolean(on_arts1_, [off_sec_, off_loc_], 'difference');

      let on_secs2_ = mfn.poly2d.Boolean(on_secs1_, [off_art_, off_loc_], 'difference');

      let on_locs2_ = mfn.poly2d.Boolean(on_locs1_, [off_art_, off_sec_], 'difference');

      let on_arts3_ = mfn.poly2d.Union(on_arts2_);

      let on_secs3_ = mfn.poly2d.Union(on_secs2_);

      let on_locs3_ = mfn.poly2d.Union(on_locs2_);

      let art_sec_ = mfn.poly2d.Boolean(off_art_, off_sec_, 'intersect');

      let art_loc_ = mfn.poly2d.Boolean(off_art_, off_loc_, 'intersect');

      let sec_loc_ = mfn.poly2d.Boolean(off_sec_, off_loc_, 'intersect');

      let art_sec1_ = mfn.poly2d.Boolean(art_sec_, block_, 'intersect');

      let art_loc1_ = mfn.poly2d.Boolean(art_loc_, block_, 'intersect');

      let sec_loc1_ = mfn.poly2d.Boolean(sec_loc_, block_, 'intersect');

      let new_parts_ = ifn.listFlat([art_sec1_, art_loc1_, sec_loc1_, on_arts3_, on_secs3_, on_locs3_, off_grids_]);

      mfn.attrib.Set(art_sec1_, `type`, 'art_sec');

      mfn.attrib.Set(art_loc1_, `type`, 'art_loc');

      mfn.attrib.Set(sec_loc1_, `type`, 'sec_loc');

      mfn.attrib.Set(on_arts3_, `type`, 'art');

      mfn.attrib.Set(on_secs3_, `type`, 'sec');

      mfn.attrib.Set(on_locs3_, `type`, 'loc');

      mfn.attrib.Set(off_grids_, `type`, 'off_grid');

      mfn.attrib.Set(new_parts_, `class`, "part");

      await exec_cluster1_node_agq2g6vtb3__transferEdgeAttribsBtwTouchingParts_($p, new_parts_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_node_agq2g6vtb3__transferEdgeAttribsTouching_($p, mfn.query.Get('_e', block_), mfn.query.Filter(mfn.query.Get('_e', new_parts_), ['road', null], '==', null));
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_node_agq2g6vtb3__copyAttribs_($p, block_, new_parts_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      mfn.edit.Delete([art_sec_, art_loc_, sec_loc_], 'delete_selected');

      mfn.edit.Delete([off_art_, off_sec_, off_loc_, on_arts1_, on_arts2_, on_secs1_, on_secs2_, on_locs1_, on_locs2_, block_], 'delete_selected');

      mfn.edit.Delete([roads_art_, roads_sec_, roads_loc_], 'delete_selected');

      return new_parts_;
    }


    async function exec_cluster1_node_agq2g6vtb3__angDot_($p, vert_) {

      let edges_ = mfn.query.Get('_e', vert_);

      let vecs_ = mfn.calc.Vector(edges_);

      let vec0_ = ifn.vecNorm(vecs_[pythonList(0, vecs_.length)]);

      let vec1_ = ifn.vecNorm(vecs_[pythonList(1, vecs_.length)]);

      return ifn.vecDot(vec0_, vec1_);
    }


    async function exec_cluster1_node_agq2g6vtb3__cleanPgons_($p, pgons_) {

      for (let pgon_ of pgons_) {

        await exec_cluster1_node_agq2g6vtb3__cleanPgonEdge_($p, pgon_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        await exec_cluster1_node_agq2g6vtb3__cleanPgonAng_($p, pgon_);
        if ($p.terminated) {
          return mfn.getModel();
        }
      }

      return mfn.query.Get('pg', pgons_);
    }


    async function exec_cluster1_node_agq2g6vtb3__cleanPgonEdge_($p, pgon_) {

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


    async function exec_cluster1_node_agq2g6vtb3__cleanPgonAng_($p, pgon_) {

      let del_posis_ = [];

      for (let vert_ of mfn.query.Get('_v', pgon_)) {

        let dot_ = await exec_cluster1_node_agq2g6vtb3__angDot_($p, vert_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (ifn.abs(dot_) > 0.9999) {

          mfn.list.Add(del_posis_, mfn.query.Get('ps', vert_), 'to_end');
        }
      }

      mfn.edit.Delete(del_posis_, 'delete_selected');
    }


    async function exec_cluster1_node_agq2g6vtb3__getSitePlines_($p, site_, road_descr_) {

      let posis_ = [];

      for (let edge_ of mfn.query.Get('_e', site_)) {

        if (mfn.attrib.Get(edge_, 'road') == road_descr_) {

          let start_end_ = mfn.query.Get('ps', edge_);

          if (ifn.len(posis_) == 0 || start_end_[pythonList(0, start_end_.length)] != posis_[pythonList(-1, posis_.length)][pythonList(-1, posis_[pythonList(-1, posis_.length)].length)]) {

            mfn.list.Add(posis_, start_end_, 'to_end');
          } else {

            mfn.list.Add(posis_[pythonList(-1, posis_.length)], start_end_[pythonList(1, start_end_.length)], 'to_end');
          }
        }
      }

      if (ifn.len(posis_) == 0) {

        return [];
      }

      if (ifn.len(posis_) > 1 && posis_[pythonList(-1, posis_.length)][pythonList(-1, posis_[pythonList(-1, posis_.length)].length)] == posis_[pythonList(0, posis_.length)][pythonList(0, posis_[pythonList(0, posis_.length)].length)]) {

        let first_list_ = ifn.listJoin(posis_[pythonList(-1, posis_.length)], posis_[pythonList(0, posis_.length)].slice(1));

        posis_[pythonList(0, posis_.length)] = first_list_;

        posis_ = posis_.slice(0, -1);
      }

      let site_plines_ = mfn.make.Polyline(posis_, 'open');

      return site_plines_;
    }


    async function exec_cluster1_node_agq2g6vtb3__transferEdgeAttribsTouching_($p, from_edges_, to_edges_) {

      for (let to_edge_ of to_edges_) {

        let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

        let from_edge_ = await exec_cluster1_node_agq2g6vtb3__touchingEdge_($p, from_edges_, cen_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (from_edge_ != null) {

          mfn.attrib.Set(to_edge_, `road`, mfn.attrib.Get(from_edge_, 'road'));
        }
      }
    }


    async function exec_cluster1_node_agq2g6vtb3__transferEdgeAttribsBtwTouchingParts_($p, parts_) {

      let edges_ = mfn.query.Filter(mfn.query.Get('_e', parts_), ['road', null], '==', null);

      for (let to_edge_ of edges_) {

        if (mfn.attrib.Get(to_edge_, 'road') == undefined) {

          let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

          let idx_ = ifn.listFind(edges_, to_edge_);

          let from_edges_ = ifn.listJoin(edges_.slice(0, idx_), edges_.slice(idx_ + 1));

          let from_edge_ = await exec_cluster1_node_agq2g6vtb3__touchingEdge_($p, from_edges_, cen_);
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


    async function exec_cluster1_node_agq2g6vtb3__copyAttribs_($p, pgon_from_, pgons_to_) {

      let block_id_ = mfn.attrib.Get(pgon_from_, 'block_id');

      let block_type_ = mfn.attrib.Get(pgon_from_, 'block_type');

      let site_ = mfn.attrib.Get(pgon_from_, 'site');

      for (let pgon_to_ of pgons_to_) {

        mfn.attrib.Set(pgon_to_, `block_id`, block_id_);

        mfn.attrib.Set(pgon_to_, `block_type`, block_type_);

        mfn.attrib.Set(pgon_to_, `site`, site_);
      }
    }


    async function exec_cluster1_node_agq2g6vtb3__touchingEdge_($p, edges_, xyz_) {

      for (let edge_ of edges_) {

        let r_ = mfn.calc.Ray(edge_);

        let d_ = ifn.distance(xyz_, ifn.project(xyz_, r_, 0));

        if (d_ < 0.01) {

          return edge_;
        }
      }

      return null;
    }

    async function exec_cluster1_node_agq2g6vtb3($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: add on grid strips art sec', '__null__')
      }


      let blocks_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', "block_corner");

      mfn.edit.Delete(blocks_, 'keep_selected');

      for (let block_ of blocks_) {

        let result_ = await exec_cluster1_node_agq2g6vtb3_processBlock_($p, block_);
        if ($p.terminated) {
          return mfn.getModel();
        }
      }
    }


    async function exec_cluster1_node_ito4lzi5zuh_createOnGridCorners_($p, block_) {

      let roads_loc_ = await exec_cluster1_node_ito4lzi5zuh__getSitePlines_($p, block_, "road_loc");
      if ($p.terminated) {
        return mfn.getModel();
      }

      let off_ = mfn.poly2d.OffsetMitre(roads_loc_, PART_LOC_D_, 100, 'square_end');

      let corner_ = mfn.poly2d.Boolean(off_, block_, 'intersect');

      let back_ = mfn.poly2d.Boolean(block_, off_, 'difference');

      mfn.attrib.Set(corner_, `type`, 'loc_loc');

      mfn.attrib.Set(corner_, `special`, 'concave_corner');

      mfn.attrib.Set(back_, `class`, 'leftover');

      mfn.attrib.Set(back_, `type`, 'concave_corner');

      let new_parts_ = ifn.listFlat([corner_, back_]);

      await exec_cluster1_node_ito4lzi5zuh__transferEdgeAttribsTouching_($p, mfn.query.Get('_e', block_), mfn.query.Get('_e', new_parts_));
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_node_ito4lzi5zuh__transferEdgeAttribsBtwTouchingParts_($p, new_parts_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await exec_cluster1_node_ito4lzi5zuh__copyAttribs_($p, block_, new_parts_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      mfn.edit.Delete([block_, off_], 'delete_selected');

      return new_parts_;
    }


    async function exec_cluster1_node_ito4lzi5zuh__angDot_($p, vert_) {

      let edges_ = mfn.query.Get('_e', vert_);

      let vecs_ = mfn.calc.Vector(edges_);

      let vec0_ = ifn.vecNorm(vecs_[pythonList(0, vecs_.length)]);

      let vec1_ = ifn.vecNorm(vecs_[pythonList(1, vecs_.length)]);

      return ifn.vecDot(vec0_, vec1_);
    }


    async function exec_cluster1_node_ito4lzi5zuh__cleanPgons_($p, pgons_) {

      for (let pgon_ of pgons_) {

        await exec_cluster1_node_ito4lzi5zuh__cleanPgonEdge_($p, pgon_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        await exec_cluster1_node_ito4lzi5zuh__cleanPgonAng_($p, pgon_);
        if ($p.terminated) {
          return mfn.getModel();
        }
      }

      return mfn.query.Get('pg', pgons_);
    }


    async function exec_cluster1_node_ito4lzi5zuh__cleanPgonEdge_($p, pgon_) {

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


    async function exec_cluster1_node_ito4lzi5zuh__cleanPgonAng_($p, pgon_) {

      let del_posis_ = [];

      for (let vert_ of mfn.query.Get('_v', pgon_)) {

        let dot_ = await exec_cluster1_node_ito4lzi5zuh__angDot_($p, vert_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (ifn.abs(dot_) > 0.9999) {

          mfn.list.Add(del_posis_, mfn.query.Get('ps', vert_), 'to_end');
        }
      }

      mfn.edit.Delete(del_posis_, 'delete_selected');
    }


    async function exec_cluster1_node_ito4lzi5zuh__getSitePlines_($p, site_, road_descr_) {

      let posis_ = [];

      for (let edge_ of mfn.query.Get('_e', site_)) {

        if (mfn.attrib.Get(edge_, 'road') == road_descr_) {

          let start_end_ = mfn.query.Get('ps', edge_);

          if (ifn.len(posis_) == 0 || start_end_[pythonList(0, start_end_.length)] != posis_[pythonList(-1, posis_.length)][pythonList(-1, posis_[pythonList(-1, posis_.length)].length)]) {

            mfn.list.Add(posis_, start_end_, 'to_end');
          } else {

            mfn.list.Add(posis_[pythonList(-1, posis_.length)], start_end_[pythonList(1, start_end_.length)], 'to_end');
          }
        }
      }

      if (ifn.len(posis_) == 0) {

        return [];
      }

      if (ifn.len(posis_) > 1 && posis_[pythonList(-1, posis_.length)][pythonList(-1, posis_[pythonList(-1, posis_.length)].length)] == posis_[pythonList(0, posis_.length)][pythonList(0, posis_[pythonList(0, posis_.length)].length)]) {

        let first_list_ = ifn.listJoin(posis_[pythonList(-1, posis_.length)], posis_[pythonList(0, posis_.length)].slice(1));

        posis_[pythonList(0, posis_.length)] = first_list_;

        posis_ = posis_.slice(0, -1);
      }

      let site_plines_ = mfn.make.Polyline(posis_, 'open');

      return site_plines_;
    }


    async function exec_cluster1_node_ito4lzi5zuh__transferEdgeAttribsTouching_($p, from_edges_, to_edges_) {

      for (let to_edge_ of to_edges_) {

        let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

        let from_edge_ = await exec_cluster1_node_ito4lzi5zuh__touchingEdge_($p, from_edges_, cen_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (from_edge_ != null) {

          mfn.attrib.Set(to_edge_, `road`, mfn.attrib.Get(from_edge_, 'road'));
        }
      }
    }


    async function exec_cluster1_node_ito4lzi5zuh__transferEdgeAttribsBtwTouchingParts_($p, parts_) {

      let edges_ = mfn.query.Filter(mfn.query.Get('_e', parts_), ['road', null], '==', null);

      for (let to_edge_ of edges_) {

        if (mfn.attrib.Get(to_edge_, 'road') == undefined) {

          let cen_ = mfn.calc.Centroid(to_edge_, 'ps_average');

          let idx_ = ifn.listFind(edges_, to_edge_);

          let from_edges_ = ifn.listJoin(edges_.slice(0, idx_), edges_.slice(idx_ + 1));

          let from_edge_ = await exec_cluster1_node_ito4lzi5zuh__touchingEdge_($p, from_edges_, cen_);
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


    async function exec_cluster1_node_ito4lzi5zuh__copyAttribs_($p, pgon_from_, pgons_to_) {

      let block_id_ = mfn.attrib.Get(pgon_from_, 'block_id');

      let block_type_ = mfn.attrib.Get(pgon_from_, 'block_type');

      let site_ = mfn.attrib.Get(pgon_from_, 'site');

      for (let pgon_to_ of pgons_to_) {

        mfn.attrib.Set(pgon_to_, `block_id`, block_id_);

        mfn.attrib.Set(pgon_to_, `block_type`, block_type_);

        mfn.attrib.Set(pgon_to_, `site`, site_);
      }
    }


    async function exec_cluster1_node_ito4lzi5zuh__touchingEdge_($p, edges_, xyz_) {

      for (let edge_ of edges_) {

        let r_ = mfn.calc.Ray(edge_);

        let d_ = ifn.distance(xyz_, ifn.project(xyz_, r_, 0));

        if (d_ < 1) {

          return edge_;
        }
      }

      return null;
    }

    async function exec_cluster1_node_ito4lzi5zuh($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: add  on grid strips concave corners', '__null__')
      }


      let blocks_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', "off_grid");

      for (let block_ of blocks_) {

        let new_parts_ = await exec_cluster1_node_ito4lzi5zuh_createOnGridCorners_($p, block_);
        if ($p.terminated) {
          return mfn.getModel();
        }

        mfn.visualize.Color(new_parts_, [1, 0, 1]);
      }
    }


    async function exec_cluster1_node_oe05a3sun4c($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: merge all cold parts', '__null__')
      }

    }


    async function exec_cluster1_node_2a54ii5s4hm($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: merge', '__null__')
      }

    }


    async function exec_cluster1_node_ieuxdb5rfmm_addAttribsCrvRoads_($p, parts_) {

      for (let part_ of parts_) {

        let part_edges_ = mfn.query.Get('_e', part_);

        let art_ = await exec_cluster1_node_ieuxdb5rfmm__getCrvRoadEdgeTypes_($p, part_edges_, "road_art");
        if ($p.terminated) {
          return mfn.getModel();
        }

        let sec_ = await exec_cluster1_node_ieuxdb5rfmm__getCrvRoadEdgeTypes_($p, part_edges_, "road_sec");
        if ($p.terminated) {
          return mfn.getModel();
        }

        let loc_ = await exec_cluster1_node_ieuxdb5rfmm__getCrvRoadEdgeTypes_($p, part_edges_, "road_loc");
        if ($p.terminated) {
          return mfn.getModel();
        }

        let road_types_ = ifn.listFlat([art_, sec_, loc_]);

        mfn.attrib.Set(part_, `class`, "part");

        let block_types_dict_ = {
          "road_art": "art",
          "road_sec": "sec",
          "road_ter": "ter",
          "road_loc": "loc"
        };

        let cats_dict_ = {
          "road_art": 1,
          "road_sec": 2,
          "road_ter": 3,
          "road_loc": 4
        };

        if (ifn.len(road_types_) == 0) {

          mfn.attrib.Set(part_, `type`, "off_grid");
        } else {
          if (ifn.len(road_types_) == 1) {

            let road_type_ = road_types_[pythonList(0, road_types_.length)];

            let block_type_ = ifn.string(block_types_dict_[pythonList(road_type_, block_types_dict_.length)]);

            mfn.attrib.Set(part_, `type`, block_type_);
          } else {

            let road_type0_ = road_types_[pythonList(0, road_types_.length)];

            let block_type0_ = ifn.string(block_types_dict_[pythonList(road_type0_, block_types_dict_.length)]);

            let road_type1_ = road_types_[pythonList(1, road_types_.length)];

            let block_type1_ = ifn.string(block_types_dict_[pythonList(road_type1_, block_types_dict_.length)]);

            mfn.attrib.Set(part_, `type`, block_type0_ + "_" + block_type1_);
          }
        }
      }
    }


    async function exec_cluster1_node_ieuxdb5rfmm__getCrvRoadEdgeTypes_($p, part_edges_, road_descr_) {

      let edges_ = [[]];

      for (let edge_ of part_edges_) {

        let edge_road_ = mfn.attrib.Get(edge_, 'road');

        if (edge_road_ == road_descr_) {

          mfn.list.Add(edges_[pythonList(-1, edges_.length)], edge_, 'to_end');
        } else {

          if (ifn.len(edges_[pythonList(-1, edges_.length)]) != 0) {

            mfn.list.Add(edges_, [], 'to_end');
          }
        }
      }

      if (ifn.len(edges_) == 2) {

        edges_ = ifn.listJoin(edges_[pythonList(1, edges_.length)], edges_[pythonList(0, edges_.length)]);
      } else {

        edges_ = edges_[pythonList(0, edges_.length)];
      }

      if (ifn.len(edges_) == 0) {

        return [];
      }

      if (ifn.len(edges_) == 1) {

        return [road_descr_];
      }

      let road_types_ = [road_descr_];

      for (let i_ of ifn.range(1, ifn.len(edges_))) {

        let ang_ = await exec_cluster1_node_ieuxdb5rfmm__vertAng_($p, mfn.query.Get('_v', edges_[pythonList(i_, edges_.length)])[pythonList(0, mfn.query.Get('_v', edges_[pythonList(i_, edges_.length)]).length)]);
        if ($p.terminated) {
          return mfn.getModel();
        }

        if (ifn.abs(ang_ - 180) > 45) {

          mfn.list.Add(road_types_, road_descr_, 'to_end');
        }
      }

      return road_types_;
    }


    async function exec_cluster1_node_ieuxdb5rfmm__vertAng_($p, vert_) {

      let edges_ = mfn.query.Get('_e', vert_);

      if (ifn.len(edges_) == 1) {

        return 0;
      }

      let vecs_ = mfn.calc.Vector(edges_);

      let vec0_ = ifn.vecRev(ifn.vecNorm(vecs_[pythonList(0, vecs_.length)]));

      let vec1_ = ifn.vecNorm(vecs_[pythonList(1, vecs_.length)]);

      return ifn.radToDeg(ifn.vecAng2(vec1_, vec0_, [0, 0, 1]));
    }

    async function exec_cluster1_node_ieuxdb5rfmm($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: merge', '__null__')
      }


      mfn.edit.Delete(mfn.query.Get('pl', null), 'delete_selected');

      mfn.attrib.Add('pg', 'string', "special");
    }


    async function exec_cluster1_node_nl92tfguq4_applyColours_($p, colors_dict_, pgons_) {

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

    async function exec_cluster1_node_nl92tfguq4($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
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

      colors_dict_bldg_["off_grid_too_small"] = [0.81, 0.9, 0.95];

      colors_dict_bldg_["too_deep"] = [0.81, 0.9, 0.95];

      colors_dict_bldg_["concave_corner"] = [0.81, 0.9, 0.95];

      colors_dict_bldg_["leftover"] = [0.81, 0.9, 0.95];

      await exec_cluster1_node_nl92tfguq4_applyColours_($p, colors_dict_bldg_, mfn.query.Get('pg', null));
      if ($p.terminated) {
        return mfn.getModel();
      }
    }


    async function exec_cluster1_node_vebs0yyerw($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: extract aux', '__null__')
      }


      let pg_ = mfn.query.Filter(mfn.query.Get('pg', null), ['inputdata', null], '==', true);

      let pl_ = mfn.query.Filter(mfn.query.Get('pl', null), ['inputdata', null], '==', true);

      let pavements_ = mfn.query.Filter(mfn.query.Get('pg', null), ['type', null], '==', 'pavement');

      let trees_pls_ = mfn.query.Filter(mfn.query.Get('pl', null), ['type', null], '==', 'tree');

      let trees_pts_ = mfn.query.Filter(mfn.query.Get('pt', null), ['type', null], '==', 'tree');

      mfn.edit.Delete([pg_, pl_, pavements_, trees_pts_], 'keep_selected');
    }


    async function exec_cluster1_node_fq25ktmf1ft_debugCheckAreas_($p, data_) {

      let check_col_l_ = await exec_cluster1_node_fq25ktmf1ft_checkTotalsColL_($p, data_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let check_col_m_ = await exec_cluster1_node_fq25ktmf1ft_checkTotalsColM_($p, data_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let check_col_o_ = await exec_cluster1_node_fq25ktmf1ft_checkTotalsColO_($p, data_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      let check_col_g_ = await exec_cluster1_node_fq25ktmf1ft_checkTotalsColG_($p, data_);
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


    async function exec_cluster1_node_fq25ktmf1ft_checkTotalsColL_($p, data_) {

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


    async function exec_cluster1_node_fq25ktmf1ft_checkTotalsColM_($p, data_) {

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


    async function exec_cluster1_node_fq25ktmf1ft_checkTotalsColO_($p, data_) {

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


    async function exec_cluster1_node_fq25ktmf1ft_checkTotalsColG_($p, data_) {

      let total_ = 0;

      total_ = total_ + data_["road_area_art"];

      total_ = total_ + data_["road_area_sec"];

      total_ = total_ + data_["road_area_loc"];

      return total_;
    }

    async function exec_cluster1_node_fq25ktmf1ft($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_) {
      if (mfn.getModel().debug) {
        printFunc($p.console, 'Executing: End', '__null__')
      }


      let data_ = {};

      data_["road_len_art_100"] = mfn.attrib.Get(null, 'road_len_art_100');

      data_["road_len_sec_100"] = mfn.attrib.Get(null, 'road_len_sec_100');

      data_["road_len_loc_100"] = mfn.attrib.Get(null, 'road_len_loc_100');

      data_["road_len_art_50"] = mfn.attrib.Get(null, 'road_len_art_50');

      data_["road_len_sec_50"] = mfn.attrib.Get(null, 'road_len_sec_50');

      data_["road_len_loc_50"] = 0;

      data_["open_art_art_area"] = 0;

      data_["open_art_sec_area"] = 0;

      data_["open_art_loc_area"] = 0;

      data_["open_art_area"] = 0;

      data_["open_sec_sec_area"] = 0;

      data_["open_sec_loc_area"] = 0;

      data_["open_sec_area"] = 0;

      data_["open_loc_loc_area"] = mfn.attrib.Get(null, 'open_loc_loc_area');

      data_["open_loc_area"] = mfn.attrib.Get(null, 'open_loc_area');

      data_["open_og_clus0_on_art_area"] = mfn.attrib.Get(null, 'open_og_clus0_on_art_area');

      data_["open_og_clus0_on_sec_area"] = mfn.attrib.Get(null, 'open_og_clus0_on_sec_area');

      data_["open_og_clus0_on_loc_area"] = mfn.attrib.Get(null, 'open_og_clus0_on_loc_area');

      data_["open_og_clus1_on_loc_area"] = mfn.attrib.Get(null, 'open_og_clus1_on_loc_area');

      data_["open_og_clus2_on_loc_area"] = mfn.attrib.Get(null, 'open_og_clus2_on_loc_area');

      data_["open_total_area"] = mfn.attrib.Get(null, 'open_total_area');

      data_["amen_art_art_area"] = 0;

      data_["amen_art_sec_area"] = 0;

      data_["amen_art_loc_area"] = 0;

      data_["amen_art_area"] = 0;

      data_["amen_sec_sec_area"] = 0;

      data_["amen_sec_loc_area"] = 0;

      data_["amen_sec_area"] = mfn.attrib.Get(null, 'amen_og_entr0_on_sec_area');
      printFunc($p.console, `data["amen_sec_area"]`, data_["amen_sec_area"]);

      data_["amen_loc_loc_area"] = mfn.attrib.Get(null, 'amen_loc_loc_area');

      data_["amen_loc_area"] = mfn.attrib.Get(null, 'amen_loc_area') + mfn.attrib.Get(null, 'amen_og_entr0_on_loc_area') + mfn.attrib.Get(null, 'amen_og_entr1_on_loc_area');

      data_["amen_og_clus0_on_art_area"] = mfn.attrib.Get(null, 'amen_og_clus0_on_art_area');

      data_["amen_og_clus0_on_sec_area"] = mfn.attrib.Get(null, 'amen_og_clus0_on_sec_area');

      data_["amen_og_clus0_on_loc_area"] = mfn.attrib.Get(null, 'amen_og_clus0_on_loc_area');

      data_["amen_og_clus1_on_loc_area"] = mfn.attrib.Get(null, 'amen_og_clus1_on_loc_area');

      data_["amen_og_clus2_on_loc_area"] = mfn.attrib.Get(null, 'amen_og_clus2_on_loc_area');

      data_["amen_total_area"] = mfn.attrib.Get(null, 'amen_total_area');

      data_["lot_art_art_area"] = mfn.attrib.Get(null, 'lot_art_art_area');

      data_["lot_art_sec_area"] = mfn.attrib.Get(null, 'lot_art_sec_area');

      data_["lot_art_loc_area"] = mfn.attrib.Get(null, 'lot_art_loc_area');

      data_["lot_art_area"] = mfn.attrib.Get(null, 'lot_art_area') + mfn.attrib.Get(null, 'og_entr0_on_art_area');

      data_["lot_sec_sec_area"] = mfn.attrib.Get(null, 'lot_sec_sec_area');

      data_["lot_sec_loc_area"] = mfn.attrib.Get(null, 'lot_sec_loc_area');

      data_["lot_sec_area"] = mfn.attrib.Get(null, 'lot_sec_area') + mfn.attrib.Get(null, 'og_entr0_on_sec_area');

      data_["lot_loc_loc_area"] = mfn.attrib.Get(null, 'lot_loc_loc_area');

      data_["lot_loc_area"] = mfn.attrib.Get(null, 'lot_loc_area') + mfn.attrib.Get(null, 'og_entr0_on_loc_area') + mfn.attrib.Get(null, 'og_entr1_on_loc_area');

      data_["og_clus0_on_art_area"] = mfn.attrib.Get(null, 'og_clus0_on_art_area');

      data_["og_clus0_on_sec_area"] = mfn.attrib.Get(null, 'og_clus0_on_sec_area');

      data_["og_clus0_on_loc_area"] = mfn.attrib.Get(null, 'og_clus0_on_loc_area');

      data_["og_clus1_on_loc_area"] = mfn.attrib.Get(null, 'og_clus1_on_loc_area');

      data_["og_clus2_on_loc_area"] = mfn.attrib.Get(null, 'og_clus2_on_loc_area');

      data_["lot_art_art_num"] = mfn.attrib.Get(null, 'lot_art_art_num');

      data_["lot_art_sec_num"] = mfn.attrib.Get(null, 'lot_art_sec_num');

      data_["lot_art_loc_num"] = mfn.attrib.Get(null, 'lot_art_loc_num');

      data_["lot_art_num"] = mfn.attrib.Get(null, 'lot_art_num');

      data_["lot_sec_sec_num"] = mfn.attrib.Get(null, 'lot_sec_sec_num');

      data_["lot_sec_loc_num"] = mfn.attrib.Get(null, 'lot_sec_loc_num');

      data_["lot_sec_num"] = mfn.attrib.Get(null, 'lot_sec_num');

      data_["lot_loc_loc_num"] = mfn.attrib.Get(null, 'lot_loc_loc_num') + mfn.attrib.Get(null, 'open_loc_loc_num') + mfn.attrib.Get(null, 'amen_loc_loc_num');

      data_["lot_loc_num"] = mfn.attrib.Get(null, 'lot_loc_num') + mfn.attrib.Get(null, 'open_loc_num') + mfn.attrib.Get(null, 'amen_loc_num');

      data_["road_area_art"] = mfn.attrib.Get(null, 'road_area_art') + mfn.attrib.Get(null, 'site_on_art_roads_area');

      data_["road_area_sec"] = mfn.attrib.Get(null, 'road_area_sec') + mfn.attrib.Get(null, 'site_on_sec_roads_area');

      data_["road_area_loc"] = mfn.attrib.Get(null, 'road_area_loc');

      data_["og_clus0_on_art_num"] = mfn.attrib.Get(null, 'og_clus0_on_art_num');

      data_["og_clus0_on_sec_num"] = mfn.attrib.Get(null, 'og_clus0_on_sec_num');

      data_["og_clus0_on_loc_num"] = mfn.attrib.Get(null, 'og_clus0_on_loc_num');

      data_["og_clus1_on_loc_num"] = mfn.attrib.Get(null, 'og_clus1_on_loc_num');

      data_["og_clus2_on_loc_num"] = mfn.attrib.Get(null, 'og_clus2_on_loc_num');

      data_["og_entr0_on_art_area"] = 0;

      data_["og_entr0_on_sec_area"] = mfn.attrib.Get(null, 'og_entr0_on_sec_area');

      data_["og_entr0_on_loc_area"] = mfn.attrib.Get(null, 'og_entr0_on_loc_area');

      data_["og_entr1_on_loc_area"] = mfn.attrib.Get(null, 'og_entr1_on_loc_area');

      data_["site_total_area"] = mfn.attrib.Get(null, 'site_total_area');

      data_["param_ogc_w"] = PART_OG_W_;

      data_["param_lot_art_d"] = PART_ART_D_;

      data_["param_lot_sec_d"] = PART_SEC_D_;

      data_["param_lot_loc_d"] = PART_LOC_D_;

      data_["param_lot_art_w"] = PLOT_ART_W_;

      data_["param_lot_sec_w"] = PLOT_SEC_W_;

      data_["param_lot_loc_w"] = PLOT_LOC_W_;

      await exec_cluster1_node_fq25ktmf1ft_debugCheckAreas_($p, data_);
      if ($p.terminated) {
        return mfn.getModel();
      }

      await mfn.io.Export(null, "cluster.sim", 'sim', 'Save to Local Storage');

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
    let ssid_exec_cluster1_node_yowzxpdo1o = mfn.model.snapshotGetActive();

    let ssid_exec_cluster1_node_6w5xaj6iu9j = ssid_exec_cluster1_node_yowzxpdo1o;

    await exec_cluster1_node_6w5xaj6iu9j($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_f9mnqqhlrz = mfn.model.snapshotNext([ssid_exec_cluster1_node_6w5xaj6iu9j]);

    await exec_cluster1_node_f9mnqqhlrz($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_8u13jv4q514 = ssid_exec_cluster1_node_f9mnqqhlrz;

    await exec_cluster1_node_8u13jv4q514($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_v300sbtwm6p = ssid_exec_cluster1_node_8u13jv4q514;

    await exec_cluster1_node_v300sbtwm6p($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_9bpb8dz1lre = mfn.model.snapshotNext([ssid_exec_cluster1_node_v300sbtwm6p]);

    await exec_cluster1_node_9bpb8dz1lre($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_1bbah65kkdi = mfn.model.snapshotNext([ssid_exec_cluster1_node_v300sbtwm6p]);

    await exec_cluster1_node_1bbah65kkdi($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_gi0wr2wsff = mfn.model.snapshotNext([ssid_exec_cluster1_node_v300sbtwm6p]);

    await exec_cluster1_node_gi0wr2wsff($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_xc8gbwkikvh = mfn.model.snapshotNext([ssid_exec_cluster1_node_6w5xaj6iu9j]);

    await exec_cluster1_node_xc8gbwkikvh($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_ytisso4wgrr = mfn.model.snapshotNext([ssid_exec_cluster1_node_xc8gbwkikvh]);

    await exec_cluster1_node_ytisso4wgrr($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_joqspe73iin = mfn.model.snapshotNext([ssid_exec_cluster1_node_ytisso4wgrr]);

    await exec_cluster1_node_joqspe73iin($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_bzwyrx9mkk6 = ssid_exec_cluster1_node_joqspe73iin;

    await exec_cluster1_node_bzwyrx9mkk6($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_q0cbraoo5qe = ssid_exec_cluster1_node_bzwyrx9mkk6;

    await exec_cluster1_node_q0cbraoo5qe($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_ua4tq4tdncr = ssid_exec_cluster1_node_q0cbraoo5qe;

    await exec_cluster1_node_ua4tq4tdncr($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_k9k87ylilf = ssid_exec_cluster1_node_ua4tq4tdncr;

    await exec_cluster1_node_k9k87ylilf($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_sqstwmuvkl8 = ssid_exec_cluster1_node_k9k87ylilf;

    await exec_cluster1_node_sqstwmuvkl8($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_cl6rfx4972q = ssid_exec_cluster1_node_sqstwmuvkl8;

    await exec_cluster1_node_cl6rfx4972q($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_zexd2z8qt7 = mfn.model.snapshotNext([ssid_exec_cluster1_node_ytisso4wgrr]);

    await exec_cluster1_node_zexd2z8qt7($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_agq2g6vtb3 = ssid_exec_cluster1_node_zexd2z8qt7;

    await exec_cluster1_node_agq2g6vtb3($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_ito4lzi5zuh = ssid_exec_cluster1_node_agq2g6vtb3;

    await exec_cluster1_node_ito4lzi5zuh($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_oe05a3sun4c = mfn.model.snapshotNext([ssid_exec_cluster1_node_cl6rfx4972q, ssid_exec_cluster1_node_ito4lzi5zuh]);

    await exec_cluster1_node_oe05a3sun4c($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_2a54ii5s4hm = mfn.model.snapshotNext([ssid_exec_cluster1_node_gi0wr2wsff, ssid_exec_cluster1_node_oe05a3sun4c]);

    await exec_cluster1_node_2a54ii5s4hm($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_ieuxdb5rfmm = mfn.model.snapshotNext([ssid_exec_cluster1_node_9bpb8dz1lre, ssid_exec_cluster1_node_1bbah65kkdi, ssid_exec_cluster1_node_2a54ii5s4hm]);

    await exec_cluster1_node_ieuxdb5rfmm($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_nl92tfguq4 = ssid_exec_cluster1_node_ieuxdb5rfmm;

    await exec_cluster1_node_nl92tfguq4($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_vebs0yyerw = mfn.model.snapshotNext([ssid_exec_cluster1_node_6w5xaj6iu9j]);

    await exec_cluster1_node_vebs0yyerw($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);

    let ssid_exec_cluster1_node_fq25ktmf1ft = mfn.model.snapshotNext([ssid_exec_cluster1_node_nl92tfguq4, ssid_exec_cluster1_node_vebs0yyerw]);

    return await exec_cluster1_node_fq25ktmf1ft($p, IN_MODEL_, FORCE_STOP_MOBIUS_, ROAD_ART_W_, ROAD_SEC_W_, ROAD_LOC_W_, PART_ART_D_, PART_SEC_D_, PART_LOC_D_, PART_OG_D_, PART_OG_W_, PLOT_ART_W_, PLOT_SEC_W_, PLOT_LOC_W_, BLK_ART_NUM_OG_D_, BLK_ART_NUM_OG_W_, BLK_SEC_NUM_OG_D_, BLK_SEC_NUM_OG_W_, BLK_LOC_NUM_OG_D_, BLK_LOC_NUM_OG_W_, PATH_W_, OPEN_PERCENT_, AMEN_PERCENT_, PAVEMENT_W_, ADD_TREES_, TREE_SPACING_, TREE_HEIGHT_START_, TREE_HEIGHT_MAX_);
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
  const result = await exec_cluster1($p, IN_MODEL, FORCE_STOP_MOBIUS, ROAD_ART_W, ROAD_SEC_W, ROAD_LOC_W, PART_ART_D, PART_SEC_D, PART_LOC_D, PART_OG_D, PART_OG_W, PLOT_ART_W, PLOT_SEC_W, PLOT_LOC_W, BLK_ART_NUM_OG_D, BLK_ART_NUM_OG_W, BLK_SEC_NUM_OG_D, BLK_SEC_NUM_OG_W, BLK_LOC_NUM_OG_D, BLK_LOC_NUM_OG_W, PATH_W, OPEN_PERCENT, AMEN_PERCENT, PAVEMENT_W, ADD_TREES, TREE_SPACING, TREE_HEIGHT_START, TREE_HEIGHT_MAX);
  if (result === mfn.getModel()) {
    return { "model": mfn.getModel(), "result": null };
  }
  return { "model": mfn.getModel(), "result": result };
  /** * **/

}

module.exports = cluster1;
