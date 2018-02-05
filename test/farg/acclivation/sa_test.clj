(ns farg.acclivation.sa-test
  (:refer-clojure :exclude [rand rand-int cond memoize])
  (:require [better-cond.core :refer [cond]]
            [clojure.edn :as edn]
            [clojure.tools.trace :refer [deftrace] :as trace]
            [clojure.pprint :refer [pprint]]
            [clojure.test :refer :all]
            [farg.acclivation.sa :refer :all]
            [farg.util :as util :refer [dd dde]]
            [farg.with-state :refer [with-state]]
            [ubergraph.core :as uber]))

(deftest test-sa
  #_(let [g (uber/digraph [:a :b -1.0])]
    (uber/pprint (spread-activation g {:a 0.5} :iterations 2))))
