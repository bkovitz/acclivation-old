(ns farg.acclivation.sa
  (:refer-clojure :exclude [rand rand-int])
  (:require [clojure.tools.trace :refer :all]
            [clojure.pprint :refer [pprint]]
            [clojure.math.numeric-tower :as math]
            [com.rpl.specter :refer :all]
            [ubergraph.core :as uber]
            [farg.util :refer [dd dde vector-contains?] :as util]
            [farg.with-state :refer [with-state]]
            [farg.x.navs :refer :all]))

(defn make-sigmoid-fn
  "Returns a logistic function with given center and range, with given
  slope at (x-center, (/ (+ y-max y-min) 2))."
  [x-center y-min y-max slope]
  (let [y-scale (- y-max y-min)
        y-center (/ (+ y-max y-min)
                    2)
        y-offset (- y-center (/ y-scale 2))]
    (fn [x]
      (+ (/ y-scale
            (+ 1.0 (Math/exp (* slope (- x-center x)))))
         y-offset))))

(def slope-for-attractor-0p5 2.1972274554893376)

(def slope-for-attractor-0p1 2.0481573722208846)
;Floating-point inaccuracy messes this one up.

(def squash
  (make-sigmoid-fn 0.0 -1.0 1.0 slope-for-attractor-0p5))
  ;Has attractors at 0.5 and -0.5.

(defn find-target-slope [target]
  (loop [slope 2.4, lo 1.0, hi 8.0]
    (let [f (make-sigmoid-fn 0.0 -1.0 1.0 slope)
          fixed-point (nth (iterate f 1.00) 100)
          diff (- fixed-point target)]
      (println slope fixed-point)
      (if (< (Math/abs diff) 0.0000000000001)
          slope
          (if (pos? diff)
              (recur (util/midpoint slope lo) lo slope)
              (recur (util/midpoint slope hi) slope hi))))))

(defn a
  "Activation level."
  [g node-or-edge]
  (uber/attr g node-or-edge :a))

(defn node-activations [g]
  (->> g (uber/nodes) (map (fn [node] [node (a g node)])) (into {})))

(defn in-nodes [g node]
  (->> (uber/in-edges g node) (map :src)))

(defn simple-output [old-activation activation]
  activation)

(defn output-delta [old-activation activation]
  (- activation old-activation))

(defn simple-newa [old-activation sum-of-inputs]
  (squash (+ old-activation sum-of-inputs)))

(defn spread-activation
  "g is an ubergraph.
  initial-activations is map of nodes to activations. Any nodes in g without
  a corresponding key in initial-activations start with activation 0.0.

  f-output takes [old-activation activation iter], returns number for all
  outgoing edges.
  f-newa takes [old-activation sum-of-inputs], returns new activation.

  Setting :watch true prints g after each iteration.

  Returns graph with updated activations, on both nodes and edges."
  [g initial-activations
   & {:keys [iterations activation-key f-output f-newa watch edge-weight-key
             decay]
      :or {iterations 1, activation-key :activation, edge-weight-key :weight,
           f-output simple-output, f-newa simple-newa, watch false,
           decay 1.0}}]
  (let [nodes (uber/nodes g)]
    (letfn [(weight [g edge]
              (uber/attr g edge edge-weight-key))
            (set-a [g node-or-edge new-a] ;set activation
              (uber/add-attr g node-or-edge activation-key new-a))
            (a [g node-or-edge] ;activation
              (uber/attr g node-or-edge activation-key))
            (apply-default-activations [g]
              (with-state [g g]
                (doseq [node nodes]
                  (when (nil? (a g node))
                        (set-a node 0.0)))))
            (apply-initial-activations [g]
              (with-state [g g]
                (doseq [[node new-a] initial-activations]
                  (set-a node new-a))))
            (print-info [g iter]
              (println "Iteration" iter)
              (uber/pprint g))
            (one-iteration [prev-g g iter]
              (with-state [g g]
                (doseq [node nodes]
                  (bind output (* (f-output (a prev-g node) (a g node))
                                  decay #_(math/expt decay iter)))
                  (doseq [edge (uber/out-edges g node)]
                    (set-a edge (* output (weight g edge)))))
                (doseq [node nodes]
                  (bind sum-of-inputs
                        (reduce + 0.0 (->> (uber/in-edges g node)
                                           (map #(a g %)))))
                  (set-a node (f-newa (a g node) sum-of-inputs)))
                (when watch -- (print-info g iter))))]
      (loop [prev-g (apply-default-activations g)
             g (apply-initial-activations prev-g)
             iter 1]
        (if (> iter iterations)
            g
            (recur g (one-iteration prev-g g iter) (inc iter)))))))

;;; The new way, as of 2-Feb-2018

(def ^{:doc "Activation"} A (attr :a))

(def decay 1.0)

(defn apply-initial-activations
  [g initial-activations]
  (with-state [g g]
    (doseq [[node new-a] initial-activations]
      (uber/add-attr node :a new-a))))

(defn incoming-to [gnode]
  (->> gnode
    (traverse [IN-EDGES (collect-one WEIGHT) SRC A])
    (reduce (fn [total [weight a]] (+ total (* weight a))) 0)))

(defn spread-activation
  [g initial-activations & {:keys [iterations] :or {iterations 1}}]
  (let [g0 (transform [NODES VAL A]
             (fn [[_ node] _]
               (get initial-activations node 0.0))
             g)]
    (reduce (fn [g iteration] (transform [NODES VAL A]
                (fn [gnode a]
                  (squash (+ a (* decay (incoming-to gnode)))))
                g))
            g0
            (range iterations))))

;IDEA Try allowing activation to feed into edge weights.
          
;  (let [g (transform [NODES VAL A] (fn [[_ node] _] (get initial-activations node 0.0)) g)]
;  #_(->> g ;(apply-initial-activations g initial-activations)
;    #_(transform [NODES VAL A] (fn [node _] (get initial-activations node 0.0)))
;    (transform [NODES VAL A]
;               (fn [gnode a]
;                 (squash (+ a (* decay (incoming-to gnode)))))))
