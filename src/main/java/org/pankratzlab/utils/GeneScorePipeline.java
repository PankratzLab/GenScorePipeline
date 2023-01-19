package org.pankratzlab.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Objects;
import java.util.Optional;
import java.util.Set;
import java.util.StringJoiner;
import java.util.Vector;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.collections4.MultiSet;
import org.apache.commons.collections4.multiset.HashMultiSet;
import org.apache.commons.math3.stat.inference.TTest;
import org.apache.commons.math3.stat.inference.WilcoxonSignedRankTest;
import org.pankratzlab.common.Aliases;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CmdLine;
import org.pankratzlab.common.Command;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.GenomeBuild;
import org.pankratzlab.common.GenomicPosition;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.Resources;
import org.pankratzlab.common.Resources.CHROMOSOME;
import org.pankratzlab.common.StrandOps;
import org.pankratzlab.common.StrandOps.AlleleOrder;
import org.pankratzlab.common.StrandOps.CONFIG;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.bioinformatics.Sequence;
import org.pankratzlab.common.filesys.Positions;
import org.pankratzlab.common.parsing.DoubleFilter;
import org.pankratzlab.common.parsing.StandardFileColumns;
import org.pankratzlab.common.plots.ForestPlot;
import org.pankratzlab.common.plots.POPULATION;
import org.pankratzlab.common.stats.Maths.COMPARISON;
import org.pankratzlab.common.stats.ProbDist;
import org.pankratzlab.common.stats.RegressionModel;
import org.pankratzlab.fileparser.AbstractColumnFilter;
import org.pankratzlab.fileparser.AbstractFileColumn;
import org.pankratzlab.fileparser.AliasedFileColumn;
import org.pankratzlab.fileparser.ColumnFilter;
import org.pankratzlab.fileparser.ColumnFilters;
import org.pankratzlab.fileparser.DataLine;
import org.pankratzlab.fileparser.DoubleWrapperColumn;
import org.pankratzlab.fileparser.FileColumn;
import org.pankratzlab.fileparser.FileParser;
import org.pankratzlab.fileparser.FileParserFactory;
import org.pankratzlab.fileparser.FixedValueColumn;
import org.pankratzlab.fileparser.IndexedFileColumn;
import org.pankratzlab.fileparser.IntegerWrapperColumn;
import org.pankratzlab.fileparser.ParseFailureException;
import org.pankratzlab.utils.MergeExtractPipeline.DataSource;
import org.pankratzlab.utils.bioinformatics.ParseSNPlocations;
import org.pankratzlab.utils.filesys.SnpMarkerSet;
import org.pankratzlab.utils.gwas.DosageData;
import org.pankratzlab.utils.gwas.DosageData.Trio;
import org.pankratzlab.utils.gwas.windows.HitWindows;

import com.google.common.base.Joiner;
import com.google.common.base.Predicates;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multiset;
import com.google.common.collect.Sets;
import com.google.common.collect.Streams;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;
import com.googlecode.charts4j.collect.Lists;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;

public class GeneScorePipeline {

  private static final String FILE_EXT_OUT = ".out";
  private static final String FILE_PREFIX_HITS = "/hits_";
  private static final String ARG_WORKING_DIR = "workDir=";
  private static final String ARG_INDEX_THRESH = "indexThresh=";
  private static final String ARG_WINDOW_SIZE = "minWinSize=";
  private static final String ARG_WINDOW_EXT = "winThresh=";
  private static final String ARG_MISS_THRESH = "missThresh=";
  private static final String ARG_R_LIBS_DIR = "rLibsDir=";
  private static final String FLAG_PLOT_ODDS_RATIO = "-plotOddsRatio";

  private static float DEFAULT_INDEX_THRESHOLD = (float) 0.00000005;
  private static int DEFAULT_WINDOW_MIN_SIZE_PER_SIDE = 500000;// 500kb each side is technically a
                                                               // 1M window until the next hit
                                                               // region, but we now take this into
                                                               // consideration in the main
                                                               // algorithm
  private static float DEFAULT_WINDOW_EXTENSION_THRESHOLD = (float) 0.000005; // (float)0.00001;
  private static String[] DEFAULT_ADDL_ANNOT_VAR_NAMES = new String[0];
  private static double DEFAULT_MIN_MISS_THRESH = 0.5;

  private static final String CROSS_FILTERED_DATAFILE = "bimData.xln";
  private static final String DATA_SOURCE_FILENAME = "data.txt";
  private static final String META_DIRECTORY = "meta_results";
  private static final String META_ANALYSIS_DEFINITIONS_FILE = "meta.analysis";

  private static final String MARKER_COL_NAME = "Marker";
  private static final String EFFECT_ALLELE_COL_NAME = "EffectAllele";
  private static final String NON_EFFECT_ALLELE_COL_NAME = "NonEffectAllele";
  private static final String BETA_COL_NAME = "beta";
  private static final String REGRESSION_HEADER = new StringJoiner("\t").add("STUDY")
                                                                        .add("DATAFILE")
                                                                        .add("INDEX-THRESHOLD")
                                                                        .add("PHENOTYPE")
                                                                        .add("SUBTYPE")
                                                                        .add("VARIATE")
                                                                        .add("BASE-R-SQR")
                                                                        .add("R-SQR").add("R-DIFF")
                                                                        .add("P-VALUE")
                                                                        .add("EXCEL-SIG")
                                                                        .add("BETA").add("SE")
                                                                        .add("DIRECTION").add("NUM")
                                                                        .add("CASES")
                                                                        .add("CONTROLS")
                                                                        .add("PAIRED-T-P-VALUE")
                                                                        .add("WILCOXON-SIGNED-RANK-P-VALUE")
                                                                        .add("PAIRED-STAT-NUM-TRIOS")
                                                                        .add("#sigInMeta")
                                                                        .add("#indexVariantsInMeta")
                                                                        .add("#indexVariantsInDataset")
                                                                        .add("B-F-SCORE")
                                                                        .add("INVCHI-SCORE")
                                                                        .toString();
  private static final String MARKER_RESULT_HEADER = new StringJoiner("\t").add("STUDY")
                                                                           .add("DATAFILE")
                                                                           .add("INDEX-THRESHOLD")
                                                                           .add("FACTOR")
                                                                           .add("MARKERNAME")
                                                                           .add("META_BETA")
                                                                           .add("META_SE")
                                                                           .add("BETA").add("SE")
                                                                           .toString();
  private final String metaDir;
  private final String rLibsDir;
  private final boolean plotOddsRatio;
  private final boolean overwrite;

  private float[] indexThresholds = new float[] {DEFAULT_INDEX_THRESHOLD};
  private int[] windowMinSizePerSides = new int[] {DEFAULT_WINDOW_MIN_SIZE_PER_SIDE};
  private float[] windowExtensionThresholds = new float[] {DEFAULT_WINDOW_EXTENSION_THRESHOLD};
  private double minMissThresh = DEFAULT_MIN_MISS_THRESH;
  private int cmacThresh = 5;

  // private int numThreads = 1;
  // private boolean runPlink = false;
  // private boolean runRegression = false;
  // private boolean writeHist = false;

  static class MetaFile {
    public String metaFile;
    public String metaRoot;
    public Map<String, MetaMarker> metaMarkers;
    public Map<String, String[]> metaMarkerAliasLookup;
    private final int hash;

    public MetaFile(String file) {
      this.metaFile = file;
      this.metaRoot = ext.rootOf(file, false);
      this.metaMarkers = new HashMap<>();
      this.metaMarkerAliasLookup = new HashMap<>();
      this.hash = Objects.hash(metaFile, metaRoot);
    }

    @Override
    public int hashCode() {
      return hash;
    }

    @Override
    public boolean equals(Object obj) {
      if (this == obj) return true;
      if (obj == null) return false;
      if (getClass() != obj.getClass()) return false;
      MetaFile other = (MetaFile) obj;
      return Objects.equals(metaFile, other.metaFile) && Objects.equals(metaRoot, other.metaRoot);
    }

  }

  static class MetaMarker {
    final String name;
    final double beta, se, pval, freq;
    final Alleles alleles;

    public MetaMarker(String name, double beta, double se, double pval, double freq, String a1,
                      String a2) {
      this.name = name;
      this.beta = beta;
      this.se = se;
      this.freq = freq;
      this.pval = pval;
      this.alleles = new Alleles(a1, a2);
    }
  }

  static class Alleles {

    final String a1;
    final String a2;

    public Alleles(String a1, String a2) {
      this.a1 = a1;
      this.a2 = a2;
    }

    @Override
    public int hashCode() {
      return Objects.hash(a1, a2);
    }

    @Override
    public boolean equals(Object obj) {
      if (this == obj) {
        return true;
      }
      if (!(obj instanceof Alleles)) {
        return false;
      }
      Alleles other = (Alleles) obj;
      return Objects.equals(a1, other.a1) && Objects.equals(a2, other.a2);
    }

  }

  private final ArrayList<MetaFile> metaFiles = new ArrayList<>();
  private final ArrayList<String> covarFiles = new ArrayList<>();
  private final ArrayList<Study> studies = new ArrayList<>();
  private final Map<String, Map<String, Double>> covarData = new HashMap<>();
  private final List<String> covarOrder = new ArrayList<>();
  private final List<Constraint> analysisConstraints = new ArrayList<>();
  private final HashMap<String, HashMap<String, Integer>> dataCounts = new HashMap<>();
  private final List<MetaAnalysis> metaAnalyses = new ArrayList<>();
  private final Table<Constraint, MetaFile, Map<MetaAnalysis, MetaAnalysisSummary>> summaryResults = HashBasedTable.create();

  // private int bimChrIndex = 0;
  // private int bimMkrIndex = 1;
  // private int bimPosIndex = 3;
  // private int bimA1Index = 4;
  // private int bimA2Index = 5;

  private final int hitsMkrIndex = 1;
  private final boolean runMetaHW;

  private final Logger log;
  private final GenomeBuild build;

  private static class HitMarker {

    private final String effectAllele;
    private final String nonEffectAllele;
    private final Double effect;

    /**
     * @param effectAllele
     * @param nonEffectAllele
     * @param effect
     */
    private HitMarker(String effectAllele, String nonEffectAllele, Double effect) {
      super();
      this.effectAllele = effectAllele;
      this.nonEffectAllele = nonEffectAllele;
      this.effect = effect;
    }

    public String getEffectAllele() {
      return effectAllele;
    }

    public String getNonEffectAllele() {
      return nonEffectAllele;
    }

    public Double getEffect() {
      return effect;
    }

  }

  private class Study {

    String studyDir;
    String studyName;

    String dataSource;
    ArrayList<DataSource> dataSources;
    // dataFile, constraint -> DosageData
    Table<Constraint, MetaFile, DosageData> data = HashBasedTable.create();
    Table<Constraint, MetaFile, Table<String, Alleles, Integer>> markerDosageIndices = HashBasedTable.create();

    ArrayList<String> phenoFiles = new ArrayList<>();

    HashMap<String, PhenoData> phenoData = new HashMap<>();

    // constraint -> datafile
    Table<Constraint, MetaFile, TrioScoreTest.Results> trioScoreTests = HashBasedTable.create();
    // constraint -> datafile -> phenofile
    Table<Constraint, MetaFile, HashMap<PhenoData, RegressionResult>> regressions = HashBasedTable.create();
    // constraint -> datafile -> phenofile -> marker
    Table<Constraint, MetaFile, Table<PhenoData, String, RegressionResult>> markerRegressions = HashBasedTable.create();
    // constraint -> datafile -> phenofile -> marker
    Table<Constraint, MetaFile, Table<PhenoData, String, RegressionResult>> multivarRegressions = HashBasedTable.create();
    // constraint -> datafile -> phenofile -> marker
    Table<Constraint, MetaFile, Table<PhenoData, String, Double>> markerMAC = HashBasedTable.create();
    // constraint -> datafile
    Table<Constraint, MetaFile, double[]> scores = HashBasedTable.create();
    // constraint -> datafile
    Table<Constraint, MetaFile, Integer> hitSnpCounts = HashBasedTable.create();
    // constraint -> datafile
    Table<Constraint, MetaFile, Integer> hitWindowCnts = HashBasedTable.create();
    // constraint -> datafile -> markerName
    Table<Constraint, MetaFile, Map<String, HitMarker>> markerData = HashBasedTable.create();
    // constraint -> datafile -> markerName
    Table<Constraint, MetaFile, Map<String, Integer>> matchedMarkerIndices = HashBasedTable.create();
    // constraint -> marker -> score[id1\tid2][matchedMarkerIndex]
    Table<Constraint, MetaFile, double[][]> matchedMarkerScores = HashBasedTable.create();
    // constraint -> marker -> dosage[id1\tid2][matchedMarkerIndex]
    // dosages are replicated here because if alleles are flipped, dosages are flipped (2 - dose)
    // so unfortunately we can't just use the raw dosage data from the data table
    Table<Constraint, MetaFile, double[][]> matchedMarkerDosages = HashBasedTable.create();
    // constraint -> MetaFile -> id1\tid2 -> score
    Table<Constraint, MetaFile, Map<String, Double>> indivMetaScores = HashBasedTable.create();
    // constraint -> id1\tid2 -> mkrRatio, 2 * cnt, cnt2
    Table<Constraint, String, double[]> indivScoreData = HashBasedTable.create();

    public boolean mkrPassesMAC(Constraint constr, MetaFile mf, PhenoData pd, String mkr) {
      double DOSAGE_THRESHOLD = cmacThresh;
      double DOSAGE_MAX_THRESHOLD = (pd.indivs.size() * 2) - cmacThresh;
      double sum = markerMAC.get(constr, mf).get(pd, mkr);
      return sum > DOSAGE_THRESHOLD && sum < DOSAGE_MAX_THRESHOLD;
    }

    public void dumpStudy() {
      Set<Constraint> allConstraints = new HashSet<>();
      allConstraints.addAll(regressions.rowKeySet());
      allConstraints.addAll(scores.rowKeySet());
      allConstraints.addAll(markerData.rowKeySet());
      allConstraints.addAll(matchedMarkerScores.rowKeySet());
      allConstraints.addAll(indivMetaScores.rowKeySet());
      allConstraints.addAll(indivScoreData.rowKeySet());

      Set<MetaFile> allMetas = new HashSet<>();
      allMetas.addAll(regressions.columnKeySet());
      allMetas.addAll(scores.columnKeySet());
      allMetas.addAll(markerData.columnKeySet());
      allMetas.addAll(indivMetaScores.columnKeySet());

      Set<String> allIndivs = new HashSet<>();
      indivMetaScores.cellSet().stream().forEach(c -> allIndivs.addAll(c.getValue().keySet()));
      allIndivs.addAll(indivScoreData.columnKeySet());

      Table<Constraint, MetaFile, Map<String, Integer>> idIndices = HashBasedTable.create();

      for (Constraint constr : allConstraints) {
        PrintWriter writer = Files.getAppropriateWriter(studyDir + constr.analysisString
                                                        + "_audit.xln");
        writer.print("FID\tIID");
        // metaFiles by SNPS (meta_SNP)
        for (MetaFile mf : allMetas) {
          writer.print("\tSCORE_" + mf.metaRoot);
        }
        for (MetaFile mf : allMetas) {
          String[][] dataIds = this.data.get(constr, mf).getIds();
          Map<String, Integer> idIndexMap = new ConcurrentHashMap<>();
          IntStream.range(0, dataIds.length).parallel()
                   .forEach(ind -> idIndexMap.put(dataIds[ind][0] + "\t" + dataIds[ind][1], ind));
          idIndices.put(constr, mf, idIndexMap);

          Map<String, HitMarker> hitMarkerData = markerData.get(constr, mf);
          for (String mkr : hitMarkerData.keySet()) {
            HitMarker hm = hitMarkerData.get(mkr);
            writer.print("\t" + mf.metaRoot + "_" + mkr
                         + (hm == null ? "" : ("_" + hm.effectAllele)));
          }
        }
        // phenoFiles by phenos/covars (pheno_pheno, pheno_covar1)
        for (String pheno : phenoData.keySet()) {
          writer.print("\t" + pheno + "_" + pheno);
          for (String covar : phenoData.get(pheno).covars) {
            writer.print("\t" + pheno + "_" + covar);
          }
        }
        writer.println();

        final int printMissingPerPheno = 5;
        Multiset<String> missingPhenos = HashMultiset.create();

        for (String indiv : allIndivs) {
          StringJoiner sj = new StringJoiner("\t");
          sj.add(indiv);
          for (MetaFile mf : allMetas) {
            sj.add(indivMetaScores.get(constr, mf).get(indiv).toString());
          }
          // metaFiles by SNPS (meta_SNP)
          for (MetaFile mf : allMetas) {
            Map<String, HitMarker> hitMarkerData = markerData.get(constr, mf);
            final double[][] ds = matchedMarkerDosages.get(constr, mf);

            for (String mkr : hitMarkerData.keySet()) {
              final Integer idInd = idIndices.get(constr, mf).get(indiv);
              final Integer mkrIndex = matchedMarkerIndices.get(constr, mf).get(mkr);
              if (idInd != null && mkrIndex != null) {
                double dosage = ds[idInd][mkrIndex];
                sj.add(Double.toString(dosage));
              } else {
                sj.add(".");
              }
            }

          }
          // phenoFiles by phenos/covars (pheno_pheno, pheno_covar1)
          for (String pheno : phenoData.keySet()) {
            PhenoData pd = phenoData.get(pheno);
            PhenoIndiv pi = pd.indivs.get(indiv);
            if (pi == null) {
              missingPhenos.add(pheno);
              if (missingPhenos.count(pheno) <= printMissingPerPheno) {
                log.reportWarning("ID " + indiv + " not found for phenotype " + pheno);
              }
              sj.add(".");
              for (@SuppressWarnings("unused")
              String covar : pd.covars) {
                sj.add(".");
              }
            } else {
              Map<String, Double> covars = pi.covars;
              sj.add(Double.toString(pi.getDepvar()));
              for (String covar : pd.covars) {
                sj.add(covars.get(covar).toString());
              }
            }
          }
          writer.println(sj.toString());
        }
        writer.close();
        missingPhenos.entrySet().stream().filter(e -> e.getCount() > printMissingPerPheno)
                     .forEachOrdered(e -> log.reportWarning(e.getCount()
                                                            + " total IDs were not found for phenotype "
                                                            + e.getElement()));

      }

      for (Constraint constr : allConstraints) {
        for (MetaFile mf : allMetas) {
          String path = getDirPath(this, mf, constr) + "EAC.xln";
          try (PrintWriter writer = Files.getAppropriateWriter(path)) {
            Table<PhenoData, String, Double> mac = markerMAC.get(constr, mf);
            List<PhenoData> pdList = new ArrayList<>(mac.rowKeySet());
            List<String> mkrList = new ArrayList<>(mac.columnKeySet());

            StringJoiner sj = new StringJoiner("\t");
            sj.add("MarkerName");
            for (PhenoData pd : pdList) {
              sj.add(pd.phenoName);
            }
            writer.println(sj.toString());

            for (String m : mkrList) {
              sj = new StringJoiner("\t");
              sj.add(m);
              for (PhenoData pd : pdList) {
                sj.add(Double.toString(mac.get(pd, m)));
              }
              writer.println(sj.toString());
            }

          }
        }
      }

    }

    /**
     * Returns a Multimap containing all markers present in the data files that are present in
     * hitMkrSet.
     *
     * @param hitMkrSet
     * @return
     */
    public Multimap<String, String[]> retrieveMarkers(MetaFile mf, Constraint constr,
                                                      Set<String> hitMkrSet) {
      Multimap<String, String[]> returnMarkers = HashMultimap.create();
      if (data.get(constr, mf).isEmpty()) {
        return returnMarkers;
      }
      String[] mkrs = data.get(constr, mf).getMarkerSet().getMarkerNames();
      String[][] alleles = data.get(constr, mf).getMarkerSet().getAlleles();
      for (int i = 0; i < mkrs.length; i++) {
        String mkr = mkrs[i];
        String[] mkrAlleles = alleles[i];
        if (hitMkrSet.contains(mkr)) {
          returnMarkers.put(mkr, mkrAlleles);
        }
      }
      return returnMarkers;
    }

    public boolean isPlinkData() {
      if (dataSources != null) {
        for (DataSource ds : dataSources) {
          if (!(ds.dataFile.endsWith(".bed") || ds.dataFile.endsWith(".ped"))) {
            return false;
          }
        }
        return true;
      }
      return false;
    }

    public void loadDataSources(MetaFile mf, Constraint constr,
                                Map<String, int[]> markerLocationMap) {
      String[] mkrsArray = markerLocationMap.keySet().toArray(new String[0]);
      SnpMarkerSet markerSet = new SnpMarkerSet(mkrsArray, false, log);
      markerSet.setBuild(build.getNCBIBuild());
      markerSet.parseSNPlocations(log);
      int[][] markerLocations = markerSet.getChrAndPositionsAsInts();
      dataSources = MergeExtractPipeline.parseDataFile(null, markerLocations, null, dataSource, 0,
                                                       log);
      if (dataSources.size() == 0) {
        // error
        log.reportError("no data sources loaded from file: " + dataSource + "; expected "
                        + markerLocations.length);
      } else {
        final String serOutput = studyDir
                                 + ext.replaceWithLinuxSafeCharacters(mf.metaRoot + "\t"
                                                                      + constr.analysisString)
                                 + "_dosage.ser";
        boolean exists = Files.exists(serOutput);
        if (exists && overwrite) {
          log.report("Serialized data file found, but overwrite flag set - reloading data from source...");
        }
        if (overwrite || !exists) {
          log.report("Loading data file " + dataSources.get(0).dataFile);
          DosageData d0 = new DosageData.DosageDataBuilder().dosageFile(dataSources.get(0).dataFile)
                                                            .dataIndexFile(dataSources.get(0).dataIndexFile)
                                                            .mapFile(dataSources.get(0).mapFile)
                                                            .idFile(dataSources.get(0).idFile)
                                                            .markers(mkrsArray)
                                                            .markerLocationMap(markerLocationMap)
                                                            .logger(log).verbose().build();
          if (dataSources.size() > 1) {
            for (int i = 1; i < dataSources.size(); i++) {
              log.report("Loading data file " + dataSources.get(i).dataFile);
              DosageData d1 = new DosageData.DosageDataBuilder().dosageFile(dataSources.get(i).dataFile)
                                                                .dataIndexFile(dataSources.get(i).dataIndexFile)
                                                                .mapFile(dataSources.get(i).mapFile)
                                                                .idFile(dataSources.get(i).idFile)
                                                                .markers(mkrsArray)
                                                                .markerLocationMap(markerLocationMap)
                                                                .logger(log).verbose().build();
              d0 = DosageData.combine(d0, d1, DosageData.COMBINE_OP.EITHER_IF_OTHER_MISSING, false,
                                      0, log);
              System.gc();
            }
          }
          if (d0.isEmpty()) {
            log.reportError("no data for key: " + mf.metaRoot + "\t" + constr.analysisString);
          }
          data.put(constr, mf, d0);
          d0.serialize(serOutput);
          System.gc();
        } else {
          log.report("" + serOutput + " already exists, loading previously loaded data for "
                     + dataSource);
          data.put(constr, mf, DosageData.load(serOutput));
        }
      }
    }
  }

  private static class RegressionResult {

    private static class Builder {

      private boolean logistic;
      private double rsq;
      private double baseRSq;
      private double pval;
      private double beta;
      private double se;
      private int num;
      private int nCases;
      private int nControls;
      private double stats;
      private int droppedSnps = 0;

      private Builder() {

      }

      private RegressionResult build() {
        return new RegressionResult(this);
      }

      /**
       * @param logistic the logistic to set
       */
      private Builder setLogistic(boolean logistic) {
        this.logistic = logistic;
        return this;
      }

      /**
       * @param rsq the rsq to set
       */
      private Builder setRsq(double rsq) {
        this.rsq = rsq;
        return this;
      }

      /**
       * @param baseRSq the baseRSq to set
       */
      private Builder setBaseRSq(double baseRSq) {
        this.baseRSq = baseRSq;
        return this;
      }

      /**
       * @param pval the pval to set
       */
      private Builder setPval(double pval) {
        this.pval = pval;
        return this;
      }

      /**
       * @param beta the beta to set
       */
      private Builder setBeta(double beta) {
        this.beta = beta;
        return this;
      }

      /**
       * @param se the se to set
       */
      private Builder setSe(double se) {
        this.se = se;
        return this;
      }

      /**
       * @param num the num to set
       */
      private Builder setNum(int num) {
        this.num = num;
        return this;
      }

      /**
       * @param nCases the nCases to set
       */
      private Builder setnCases(int nCases) {
        this.nCases = nCases;
        return this;
      }

      /**
       * @param nControls the nControls to set
       */
      private Builder setnControls(int nControls) {
        this.nControls = nControls;
        return this;
      }

      /**
       * @param stats the stats to set
       */
      private Builder setStats(double stats) {
        this.stats = stats;
        return this;
      }

      private Builder setDroppedSnps(int num) {
        this.droppedSnps = num;
        return this;
      }

    }

    private final boolean logistic;
    private final double rsq;
    private final double baseRSq;
    private final double pval;
    private final double beta;
    private final double se;
    private final int num;
    private final int nCases;
    private final int nControls;
    private final double stats;
    private final int droppedSnps;

    private RegressionResult(Builder builder) {
      super();
      logistic = builder.logistic;
      rsq = builder.rsq;
      baseRSq = builder.baseRSq;
      pval = builder.pval;
      beta = builder.beta;
      se = builder.se;
      num = builder.num;
      nCases = builder.nCases;
      nControls = builder.nControls;
      stats = builder.stats;
      droppedSnps = builder.droppedSnps;
    }

    /**
     * @return the logistic
     */
    boolean isLogistic() {
      return logistic;
    }

    /**
     * @return the rsq
     */
    double getRsq() {
      return rsq;
    }

    /**
     * @return the baseRSq
     */
    double getBaseRSq() {
      return baseRSq;
    }

    /**
     * @return the pval
     */
    double getPval() {
      return pval;
    }

    /**
     * @return the beta
     */
    double getBeta() {
      return beta;
    }

    /**
     * @return the se
     */
    double getSe() {
      return se;
    }

    /**
     * @return the num
     */
    int getNum() {
      return num;
    }

    /**
     * @return the nCases
     */
    int getnCases() {
      return nCases;
    }

    /**
     * @return the nControls
     */
    int getnControls() {
      return nControls;
    }

    /**
     * @return the stats
     */
    double getStats() {
      return stats;
    }

    private static Builder builder() {
      return new Builder();
    }

    private static RegressionResult dummy() {
      return builder().setLogistic(true).setRsq(Double.NaN).setBaseRSq(Double.NaN)
                      .setPval(Double.NaN).setBeta(Double.NaN).setSe(Double.NaN).setNum(0)
                      .setnCases(0).setnControls(0).setStats(Double.NaN).build();
    }

  }

  private static class TrioScoreTest {

    private static class Results {

      private final int trioCount;
      private final double wilcoxonPVal;
      private final double pairedTPVal;

      private Results(TrioScoreTest trioScoreTest) {
        super();
        double[] caseScoresArray = Doubles.toArray(trioScoreTest.caseScores);
        double[] parentMeansArray = Doubles.toArray(trioScoreTest.parentScores);
        if (caseScoresArray.length == parentMeansArray.length) {
          trioCount = caseScoresArray.length;
        } else {
          throw new IllegalStateException("Case and mean parent score arrays must be parallel");
        }
        this.wilcoxonPVal = new WilcoxonSignedRankTest().wilcoxonSignedRankTest(caseScoresArray,
                                                                                parentMeansArray,
                                                                                false);
        this.pairedTPVal = new TTest().pairedTTest(caseScoresArray, parentMeansArray);
      }

      /**
       * @return the number of trios used in the paired stats calculated
       */
      public int getTrioCount() {
        return trioCount;
      }

    }

    private final List<Double> caseScores;
    private final List<Double> parentScores;

    public TrioScoreTest() {
      caseScores = new ArrayList<>();
      parentScores = new ArrayList<>();
    }

    public void addScore(double caseScore, double fatherScore, double motherScore) {
      caseScores.add(caseScore);
      parentScores.add((fatherScore + motherScore) / 2.0);
    }

    public Results runTests() {
      if (caseScores.size() < 2) return null;
      return new Results(this);
    }

  }

  private static class Constraint {

    final float indexThreshold;
    final int windowMinSizePerSide;
    final float windowExtensionThreshold;
    final String analysisString;
    final int hash;

    public Constraint(float i, int m, float w) {
      indexThreshold = i;
      windowMinSizePerSide = m;
      windowExtensionThreshold = w;
      analysisString = new StringBuilder().append(ext.formSciNot(indexThreshold, 4, false))
                                          .append("_")
                                          .append(ext.formSciNot(windowMinSizePerSide, 4, false))
                                          .append("_")
                                          .append(ext.formSciNot(windowExtensionThreshold, 4,
                                                                 false))
                                          .toString();
      hash = Objects.hash(indexThreshold, windowExtensionThreshold, windowMinSizePerSide);
    }

    @Override
    public int hashCode() {
      return hash;
    }

    @Override
    public boolean equals(Object obj) {
      if (this == obj) return true;
      if (obj == null) return false;
      if (getClass() != obj.getClass()) return false;
      Constraint other = (Constraint) obj;
      return Float.floatToIntBits(indexThreshold) == Float.floatToIntBits(other.indexThreshold)
             && Float.floatToIntBits(windowExtensionThreshold) == Float.floatToIntBits(other.windowExtensionThreshold)
             && windowMinSizePerSide == other.windowMinSizePerSide;
    }

  }

  private static class PhenoData {

    String phenoName;
    HashMap<String, PhenoIndiv> indivs = new HashMap<>();
    ArrayList<String> covars = new ArrayList<>();

    @Override
    public int hashCode() {
      return Objects.hash(phenoName);
    }

    @Override
    public boolean equals(Object obj) {
      if (this == obj) {
        return true;
      }
      if (!(obj instanceof PhenoData)) {
        return false;
      }
      PhenoData other = (PhenoData) obj;
      return Objects.equals(phenoName, other.phenoName);
    }

  }

  private class PhenoIndiv {

    private final String fid;
    private final String iid;
    private final double depvar;
    private final Map<String, Double> covars;

    /**
     * @param fid
     * @param iid
     * @param depvar
     */
    public PhenoIndiv(String fid, String iid, double depvar) {
      super();
      this.fid = fid;
      this.iid = iid;
      this.depvar = depvar;
      this.covars = new HashMap<>();
    }

    /**
     * @return the fid
     */
    public String getFid() {
      return fid;
    }

    /**
     * @return the iid
     */
    public String getIid() {
      return iid;
    }

    /**
     * @return the depvar
     */
    public double getDepvar() {
      return depvar;
    }

    /**
     * @return an unmodifiable view of the covars
     */
    public Map<String, Double> getCovars() {
      return Collections.unmodifiableMap(covars);
    }

    public void addCovar(String name, double value) {
      if (covars.putIfAbsent(name, value) != null) {
        log.reportWarning("Duplicate covars named '" + name + "', only one will be used");
      }
    }

  }

  private static HashMap<String, Double> get1000GFreq(HashMap<String, int[]> markerMap,
                                                      POPULATION population, GenomeBuild build,
                                                      Logger log) {
    POPULATION pop = population == null ? POPULATION.ALL : population;
    HashMap<Integer, HashSet<String>> mkrsByChr, nonRSByChr;
    HashSet<String> chrMkrs, nonRS;
    HashMap<String, Double> markerFreqs = new HashMap<>();

    mkrsByChr = new HashMap<>();
    nonRSByChr = new HashMap<>();
    for (java.util.Map.Entry<String, int[]> entry : markerMap.entrySet()) {
      chrMkrs = mkrsByChr.get(entry.getValue()[0]);
      nonRS = nonRSByChr.get(entry.getValue()[0]);
      if (chrMkrs == null) {
        chrMkrs = new HashSet<>();
        mkrsByChr.put(entry.getValue()[0], chrMkrs);
      }
      if (nonRS == null) {
        nonRS = new HashSet<>();
        nonRSByChr.put(entry.getValue()[0], nonRS);
      }
      chrMkrs.add(entry.getKey());
      if (!entry.getKey().startsWith("rs")) {
        nonRS.add(entry.getKey());
      }
    }

    for (Integer chr : mkrsByChr.keySet()) {
      // skip this chromosome if we don't have any markers for which we need freqencies
      if (mkrsByChr.get(chr) == null || mkrsByChr.get(chr).isEmpty()) {
        continue;
      }
      log.report("Retrieving frequencies for " + mkrsByChr.get(chr).size()
                 + " markers on chromosome " + chr);

      String afFile = Resources.genome(build, log)
                               .chr(CHROMOSOME.valueOf("C" + Integer.toString(chr)))
                               .getG1Kphase3v5AlleleFreq().get();

      FileColumn<String> snpCol = new AliasedFileColumn("SNP", "ID");
      FileColumn<Byte> chrCol = StandardFileColumns.chr("CHROM");
      FileColumn<Integer> posCol = StandardFileColumns.pos("POS");
      FileColumn<Double> afCol = DoubleWrapperColumn.wrap(new AliasedFileColumn(pop.name(),
                                                                                pop.colName));

      chrMkrs = mkrsByChr.get(chr);
      nonRS = nonRSByChr.get(chr);
      long parseFails = 0;

      try (FileParser parser = FileParserFactory.setup(afFile, snpCol, posCol, chrCol, afCol)
                                                .build()) {
        for (DataLine line : parser) {
          String snp = line.getString(snpCol);
          try {
            if (snp.startsWith("rs")) {
              if (chrMkrs.contains(snp)) {
                markerFreqs.put(snp, line.get(afCol));
                chrMkrs.remove(snp);
              }
            } else if (nonRS.contains(snp)) {
              markerFreqs.put(snp, line.get(afCol));
              nonRS.remove(snp);
            }
          } catch (ParseFailureException e) {
            parseFails++;
          }
        }
      } catch (IOException e) {
        e.printStackTrace();
      }

      if (chrMkrs.size() > 0 || nonRS.size() > 0 || parseFails > 0) {
        log.report("Chromosome " + chr + " -- Couldn't find " + chrMkrs.size() + " RS-ids and "
                   + nonRS.size() + " non-RS snps."
                   + (parseFails > 0 ? "  Failed to parse allele frequencies for " + parseFails
                                       + " markers."
                                     : ""));
      }
    }

    return markerFreqs;
  }

  // public static void preprocessLiftOver(String file, GenomeBuild fromBuild, GenomeBuild toBuild,
  // Logger log) {
  // String chainFile = Resources.genome(fromBuild, log).getLiftOverChainFile(toBuild).get();
  // if (chainFile == null) {
  // log.reportError("No chain file available for " + fromBuild.getUCSCVersion() + " to "
  // + toBuild.getUCSCVersion() + "; aborting...");
  // return;
  // }
  // LiftOver lift = new LiftOver(new File(chainFile));
  //
  // float indexThreshold;
  // int windowMinSizePerSide;
  // float windowExtensionThreshold;
  // String[][] results = HitWindows.determine(file, indexThreshold, windowMinSizePerSide,
  // windowExtensionThreshold,
  // DEFAULT_ADDL_ANNOT_VAR_NAMES, new Logger());
  //
  // if (results == null) {
  // log.reportError("HitWindows result was null for " + file + ".");
  // } else {
  // log.report("Found " + results.length + " hit windows");
  // for (int i = 1; i < results.length; i++) {
  // String[] result = results[i];
  // try {
  // hitMkrLocations.put(result[1],
  // new int[] {Integer.parseInt(result[2]), Integer.parseInt(result[3])});
  // } catch (NumberFormatException nfe) {
  // log.reportError("Failed to parse position for result " + ArrayUtils.toStr(result));
  // hitMkrLocations.put(result[1], new int[] {-1, -1});
  // }
  // }
  // }
  //
  // }

  public static Map<String, int[]> processSNPsToPositions(String[] snps, GenomeBuild build,
                                                          Logger log) {
    Set<String> parseable = new HashSet<>();
    Set<String> lookup = new HashSet<>();

    Map<String, int[]> positions = new HashMap<>();

    for (String snp : snps) {
      if (snp.contains(":")) {
        parseable.add(snp);
      } else /* if (snp.startsWith("rs")) */ {
        lookup.add(snp);
      }
    }

    for (String snp : parseable) {
      String[] pts = snp.split(":");
      if (pts.length == 2 && pts[1].contains("_")) {
        String[] pts1 = pts[1].split("_");
        pts[1] = pts1[0];
      }
      int pos = -1;
      try {
        pos = Integer.parseInt(pts[1]);
      } catch (NumberFormatException e) {
        try {
          pos = (int) Double.parseDouble(pts[1]);
        } catch (NumberFormatException e1) {
          e.printStackTrace();
        }
      }
      positions.put(snp, new int[] {Positions.chromosomeNumber(pts[0]), pos});
    }

    if (lookup.size() > 0) {
      String[] snps1 = lookup.toArray(new String[0]);
      int[][] results = ParseSNPlocations.parseSNPLocations(snps1, build, log);
      for (int i = 0; i < snps1.length; i++) {
        positions.put(snps1[i], new int[] {results[i][0], results[i][1]});
      }
    }

    return positions;
  }

  public static void preprocessDataFiles(String[] files, POPULATION pop, GenomeBuild build,
                                         GenomeBuild liftToBuild, Logger log) {
    BufferedReader reader;
    String temp, delimiter;
    String[] header, snps = null, line;
    String[][] factors;
    int[] indices;
    HashMap<String, int[]> markerMap;
    Hashtable<String, Vector<String>> fileData;
    HashMap<String, Double> freqs;

    factors = new String[][] {Aliases.MARKER_NAMES, Aliases.CHRS, Aliases.POSITIONS,
                              Aliases.PVALUES, Aliases.ALLELE_FREQS, Aliases.EFFECTS};
    for (String filename : files) {
      if (filename.endsWith(".meta")) {
        log.reportError("Cannot process a '.meta' file (" + filename
                        + "); please rename file extension.");
        continue;
      }
      try {
        reader = Files.getAppropriateReader(filename);
        temp = reader.readLine();
        delimiter = ext.determineDelimiter(temp);
        header = temp.trim().split(delimiter);
        indices = ext.indexFactors(factors, header, false, false, true, true, log);
        markerMap = new HashMap<>();
        String errorMsg = "";
        if (indices[0] == -1) {
          errorMsg = "ERROR - no MarkerName column found";
          // ERROR - couldn't find MarkerName column! COMPLETE FAIL
        }
        if (indices[5] == -1) {
          // NO BETAS! COMPLETE FAIL
          errorMsg = errorMsg.equals("") ? "ERROR - no Beta/Effect column found"
                                         : errorMsg + "; no Beta/Effect column found";
        }
        if (indices[3] == -1) {
          // NO PVALUES! COMPLETE FAIL
          errorMsg = errorMsg.equals("") ? "ERROR - no P-Value column found"
                                         : errorMsg + "; no P-Value column found";
        }
        if (errorMsg.equals("")) {
          if (indices[1] == -1 || indices[2] == -1) {
            log.report("Chromosome and/or position are missing; parsing from lookup ... ");
            // No chromosomes/positions
            snps = HashVec.loadFileToStringArray(filename, true, new int[] {indices[0]}, false);
            Map<String, int[]> positions = processSNPsToPositions(snps, build, log);
            markerMap.putAll(positions);
          } else {
            fileData = HashVec.loadFileToHashVec(filename, indices[0],
                                                 new int[] {indices[1], indices[2]}, "\t", true,
                                                 false);
            for (String key : fileData.keySet()) {
              markerMap.put(key,
                            new int[] {Positions.chromosomeNumber(fileData.get(key).get(0)
                                                                          .split("\t")[0]),
                                       Integer.parseInt(fileData.get(key).get(0).split("\t")[1])});
            }
            snps = fileData.keySet().toArray(new String[fileData.size()]);
            fileData = null;
          }
          if (indices[4] == -1) {
            // no frequencies
            log.report("Allele frequencies are missing; parsing from 1000G freqs ... ");
            freqs = get1000GFreq(markerMap, pop, build, log);
            // freqs = null;
          } else {
            freqs = null;// new HashMap<String, Double>();
            // while ((temp = reader.readLine()) != null) {
            // line = temp.trim().split(PSF.Regex.GREEDY_WHITESPACE);
            // freqs.put(line[indices[0]], Double.valueOf(line[indices[4]]));
            // }
          }
          if (build != liftToBuild) {
            log.report("LiftOver requested from " + build.getUCSCVersion() + " to "
                       + liftToBuild.getUCSCVersion());
            String chainFile = Resources.genome(build, log).getLiftOverChainFile(liftToBuild).get();
            if (chainFile == null || !Files.exists(chainFile)) {
              log.reportError("Chain file for " + build.getUCSCVersion() + " to "
                              + liftToBuild.getUCSCVersion()
                              + " not found; unable to perform LiftOver.");
            } else {
              LiftOver liftOver = new LiftOver(new File(chainFile));
              Set<String> snpSet = markerMap.keySet();
              Set<String> failed = new HashSet<>();
              long liftCount = 0;
              for (String snp : snpSet) {
                int[] chrPos = markerMap.get(snp);
                Interval lifted = liftOver.liftOver(new Interval("chr"
                                                                 + Positions.chromosomeNumberInverse(chrPos[0]),
                                                                 chrPos[1], chrPos[1]));
                if (lifted != null && lifted.getStart() != chrPos[1]) {
                  markerMap.put(snp, new int[] {Positions.chromosomeNumber(lifted.getContig()),
                                                lifted.getStart()});
                  liftCount++;
                } else if (lifted == null) {
                  failed.add(snp);
                }
              }
              log.report("Completed LiftOver on " + liftCount + " SNPs"
                         + (failed.size() > 0 ? ("; failed to lift " + failed.size() + " SNPs"
                                                 + (failed.size() < 10 ? " ("
                                                                         + ArrayUtils.toStr(failed,
                                                                                            ",")
                                                                         + ")"
                                                                       : ""))
                                              : ""));
            }
          }

          StringBuilder newHeaderSB = new StringBuilder("SNP\tChr\tPos\tFreq\tP\tBeta");
          for (int i = 0; i < header.length; i++) {
            if (i != indices[0] && i != indices[1] && i != indices[2] && i != indices[3]
                && i != indices[4] && i != indices[5]) {
              newHeaderSB.append("\t").append(header[i]);
            }
          }

          PrintWriter metaWriter = Files.getAppropriateWriter(ext.rootOf(filename, false)
                                                              + ".meta");
          metaWriter.println(newHeaderSB.toString());
          while ((temp = reader.readLine()) != null) {
            if (temp.isEmpty()) continue;
            line = temp.trim().split(PSF.Regex.GREEDY_WHITESPACE);

            String snp = line[indices[0]];
            String chr = (indices[1] == -1 || build != liftToBuild) ? "" + markerMap.get(snp)[0]
                                                                    : line[indices[1]];
            String pos = (indices[2] == -1 || build != liftToBuild) ? "" + markerMap.get(snp)[1]
                                                                    : line[indices[2]];
            String pval = line[indices[3]];
            String freq = indices[4] == -1 ? "" + (freqs == null || freqs.isEmpty()
                                                   || freqs.get(snp) == null ? 0.0 : freqs.get(snp))
                                           : line[indices[4]];
            String beta = line[indices[5]];
            StringBuilder writeLineSB = new StringBuilder();
            writeLineSB.append(snp).append("\t").append(chr).append("\t").append(pos).append("\t")
                       .append(freq).append("\t").append(pval).append("\t").append(beta);
            for (int i = 0; i < line.length; i++) {
              if (i != indices[0] && i != indices[1] && i != indices[2] && i != indices[3]
                  && i != indices[4] && i != indices[5]) {
                writeLineSB.append("\t").append(line[i]);
              }
            }
            metaWriter.println(writeLineSB.toString());
          }
          metaWriter.flush();
          metaWriter.close();
        } else {
          log.reportError(errorMsg);
          reader.close();
          continue;
        }
      } catch (IOException e) {
        e.printStackTrace();
      }
    }
  }

  public GeneScorePipeline(String workDir, float[] indexThresholds, int[] windowMins,
                           float[] windowExtThresholds, double missThresh, int cmacThresh,
                           boolean runMetaHW, String rLibsDir, boolean plotOddsRatio,
                           GenomeBuild build, boolean overwrite, Logger log) {
    this.log = log;
    this.overwrite = overwrite;
    this.build = build;
    this.metaDir = ext.verifyDirFormat(workDir);
    this.runMetaHW = runMetaHW;
    this.rLibsDir = rLibsDir;
    this.plotOddsRatio = plotOddsRatio;
    this.minMissThresh = missThresh;
    this.cmacThresh = cmacThresh;
    this.indexThresholds = indexThresholds;
    windowMinSizePerSides = windowMins;
    windowExtensionThresholds = windowExtThresholds;
    if (indexThresholds.length != windowExtThresholds.length) {
      throw new IllegalArgumentException();
    }
    setFilePrefices();
    loadStudyFolders();
    // instantiate inner hashmaps:
    for (Study study : studies) {
      for (Constraint constr : analysisConstraints) {
        for (MetaFile mf : metaFiles) {
          study.regressions.put(constr, mf, new HashMap<>());
          study.markerRegressions.put(constr, mf, HashBasedTable.create());
          study.multivarRegressions.put(constr, mf, HashBasedTable.create());
          study.markerMAC.put(constr, mf, HashBasedTable.create());
          study.indivMetaScores.put(constr, mf, new HashMap<>());
        }
      }
    }
    loadDataCounts();
    loadCovarData();
    loadMetaFileData();
    loadMetaAnalysisFile();
    runMetaHitWindowsAndLoadData();
  }

  private void setFilePrefices() {
    for (float i : indexThresholds) {
      for (int m : windowMinSizePerSides) {
        for (float w : windowExtensionThresholds) {
          final float indexThresh = i;
          final int minSize = m;
          final float windowThresh;
          if (w < indexThresh) {
            log.reportError("Window extension threshold (" + w
                            + ") is more stringent than index threshold (" + indexThresh
                            + "), index threshold will be used as window extension threshold");
            windowThresh = indexThresh;
          } else {
            windowThresh = w;
          }
          analysisConstraints.add(new Constraint(indexThresh, minSize, windowThresh));
        }
      }
    }
  }

  private void loadStudyFolders() {
    File dir = new File(metaDir);
    File[] fs = dir.listFiles();
    for (File f : fs) {
      if (f.isDirectory()) {
        // ignore meta results directory
        if (f.getName().equals(META_DIRECTORY)) continue;
        Study study = new Study();
        study.studyName = ext.rootOf(f.getAbsolutePath(), true);
        study.studyDir = f.getAbsolutePath() + "/";
        for (File f1 : f.listFiles()) {
          if (f1.getName().endsWith(".pheno")) {
            study.phenoFiles.add(f1.getName());
          }
          if (f1.getName().equals(DATA_SOURCE_FILENAME)) {
            study.dataSource = f1.getAbsolutePath();
          }
        }
        if (study.dataSource == null) {
          log.reportError("data source file {" + DATA_SOURCE_FILENAME + "} missing for study "
                          + study.studyName);
        }
        studies.add(study);
      } else if (f.getAbsolutePath().endsWith(".meta")) {
        metaFiles.add(new MetaFile(f.getName()));
      } else if (f.getAbsolutePath().endsWith(".covar")) {
        covarFiles.add(f.getName());
      }
    }
    metaFiles.sort(new Comparator<MetaFile>() {
      @Override
      public int compare(MetaFile o1, MetaFile o2) {
        return ext.removeDirectoryInfo(o1.metaRoot).compareTo(ext.removeDirectoryInfo(o2.metaRoot));
      }
    });
    covarFiles.sort(new Comparator<String>() {
      @Override
      public int compare(String o1, String o2) {
        return ext.rootOf(o1, true).compareTo(ext.rootOf(o2, true));
      }
    });
  }

  private void loadDataCounts() {
    String countsFile = metaDir + "data.cnt";

    boolean exists = (new File(countsFile).exists());

    if (exists && overwrite) {
      log.report("Data counts exists but overwrite flag set - reloading data counts...");
    }
    if (!overwrite && exists) {
      try {
        BufferedReader reader = Files.getAppropriateReader(countsFile);
        String line = null;
        while ((line = reader.readLine()) != null && !"".equals(line)) {
          String[] temp = line.split("\t");
          HashMap<String, Integer> dFileCnts = dataCounts.get(temp[0]);
          if (dFileCnts == null) {
            dFileCnts = new HashMap<>();
            dataCounts.put(temp[0], dFileCnts);
          }
          dFileCnts.put(temp[1], Integer.parseInt(temp[2]));
        }
      } catch (NumberFormatException e) {
        log.reportException(e);
        e.printStackTrace();
      } catch (IOException e) {
        log.reportIOException(countsFile);
        log.reportException(e);
        e.printStackTrace();
      }
    } else {

      for (MetaFile mf : metaFiles) {
        String dFile = mf.metaFile;
        HashSet<Constraint> threshNeed = new HashSet<>();

        HashMap<String, Integer> cnts = dataCounts.get(dFile);
        if (cnts == null) {
          threshNeed.addAll(analysisConstraints);
          cnts = new HashMap<>();
          dataCounts.put(dFile, cnts);
        } else {
          for (Constraint pref : analysisConstraints) {
            if (!cnts.containsKey(pref.analysisString)) {
              threshNeed.add(pref);
            }
          }
        }

        if (threshNeed.size() > 0) {
          try {
            BufferedReader reader = Files.getAppropriateReader(metaDir + dFile);
            String line = reader.readLine();
            String[] dataHdrs = line.split(PSF.Regex.GREEDY_WHITESPACE);
            int[] indices = ext.indexFactors(Aliases.PVALUES, dataHdrs, false, log, false);
            int ind = -1;
            for (int i : indices) {
              if (i > 0) {
                ind = i;
                break;
              }
            }
            if (ind > 0) {
              while ((line = reader.readLine()) != null) {
                if (line.isEmpty()) continue;
                String[] parts = line.split("\t");
                if (!ext.isMissingValue(parts[ind]) && ext.isValidDouble(parts[ind])) {
                  double pval = Double.parseDouble(parts[ind]);
                  for (Constraint constraint : threshNeed) {
                    if (pval < constraint.indexThreshold) {
                      Integer cnt = dataCounts.get(dFile).get(constraint.analysisString);
                      if (cnt == null) {
                        cnt = 0;
                      }
                      cnt = cnt + 1;
                      dataCounts.get(dFile).put(constraint.analysisString, cnt);
                    }
                  }
                }
              }
            }
          } catch (NumberFormatException e) {
            log.reportException(e);
            e.printStackTrace();
          } catch (IOException e) {
            log.reportIOException(metaDir + dFile);
            log.reportException(e);
            e.printStackTrace();
          }
        }
      }

      StringBuilder output = new StringBuilder();
      for (java.util.Map.Entry<String, HashMap<String, Integer>> entry : dataCounts.entrySet()) {
        for (java.util.Map.Entry<String, Integer> subEntry : entry.getValue().entrySet()) {
          output.append(entry.getKey()).append("\t").append(subEntry.getKey()).append("\t")
                .append(subEntry.getValue()).append("\n");
        }
      }
      Files.write(output.toString(), countsFile);
    }
  }

  private void loadCovarData() {
    for (String cFile : covarFiles) {
      String[] hdr = Files.getHeaderOfFile(metaDir + cFile, log);
      if (hdr.length == 2) {
        // TODO throw error
        continue;
      }
      FileColumn<?>[] data = new FileColumn[hdr.length];

      for (int i = 0; i < hdr.length; i++) {
        if (i < 2) {
          data[i] = new AliasedFileColumn(hdr[i]);
          continue;
        }
        data[i] = DoubleWrapperColumn.wrap(new AliasedFileColumn(hdr[i]));
        if (covarData.containsKey(hdr[i])) {
          String orig = hdr[i];
          int cnt = 1;
          while (covarData.containsKey(hdr[i])) {
            hdr[i] = orig + "_" + cnt;
            cnt++;
          }
          log.reportWarning("Duplicate covariate name found {" + orig
                            + "}.  Covar will be renamed to " + hdr[i]);
        }
        Map<String, Double> vars = new HashMap<>();
        covarData.put(hdr[i], vars);
        covarOrder.add(hdr[i]);
      }
      try (FileParser parser = FileParserFactory.setup(metaDir + cFile, data).build()) {
        for (DataLine line : parser) {
          String fid = line.getUnsafe((AliasedFileColumn) data[0]);
          String iid = line.getUnsafe((AliasedFileColumn) data[1]);
          for (int i = 2; i < data.length; i++) {
            covarData.get(hdr[i])
                     .put(fid + "\t" + iid,
                          line.hasValid(data[i]) ? line.getUnsafe((DoubleWrapperColumn) data[i])
                                                 : Double.NaN);
          }
        }
      } catch (IOException e) {
        log.reportIOException(metaDir + cFile);
        log.reportException(e);
        e.printStackTrace();
      }

    }
  }

  private void loadMetaFileData() {
    String[][] factors = new String[][] {Aliases.MARKER_NAMES, Aliases.EFFECTS,
                                         Aliases.ALLELE_FREQS, Aliases.PVALUES, Aliases.STD_ERRS,
                                         Aliases.ALLELES[0], Aliases.ALLELES[1]};
    for (MetaFile mf : metaFiles) {
      String metaFile = mf.metaFile;
      String fullPath = metaDir + metaFile;
      try {
        BufferedReader reader = Files.getAppropriateReader(fullPath);
        String line = reader.readLine();
        String[] temp = line.split(PSF.Regex.GREEDY_WHITESPACE, -1);
        int[] indices = ext.indexFactors(factors, temp, false, false, true, true, new Logger());
        int index = 1;
        while ((line = reader.readLine()) != null) {
          index++;
          if (line.isEmpty()) continue;
          temp = line.split(PSF.Regex.GREEDY_WHITESPACE);
          if (temp.length <= indices[1] || temp.length <= indices[2] || temp.length <= indices[3]) {
            throw new IllegalStateException("too few columns on line " + index + " of meta file "
                                            + metaFile);
          }
          if ((indices[1] != -1 && ext.isMissingValue(temp[indices[1]]))
              || ext.isMissingValue(temp[indices[2]]) || ext.isMissingValue(temp[indices[3]])) {
            continue;
          }
          String mkr = temp[indices[0]];
          double beta = indices[1] == -1 ? Double.NaN : Double.parseDouble(temp[indices[1]]);
          double freq = Double.parseDouble(temp[indices[2]]);
          double pval = Double.parseDouble(temp[indices[3]]);
          double se = indices[4] == -1 ? Double.NaN : Double.parseDouble(temp[indices[4]]);
          String a1 = indices[5] == -1 ? "." : temp[indices[5]];
          String a2 = indices[6] == -1 ? "." : temp[indices[6]];
          mf.metaMarkers.put(mkr, new MetaMarker(mkr, beta, se, pval, freq, a1, a2));
        }
      } catch (IOException e) {
        log.reportIOException(fullPath);
        log.reportException(e);
        e.printStackTrace();
      }
    }

  }

  private static class MetaAnalysisEntry {
    Study study;
    String phenoFile;

    public MetaAnalysisEntry(Study study, String phenoFile) {
      this.study = study;
      this.phenoFile = phenoFile;
    }

  }

  private static class MetaAnalysis {
    String name;
    List<MetaAnalysisEntry> entries;

    public MetaAnalysis() {
      this.entries = new ArrayList<>();
    }

    public void addEntry(String name, Study stud, String phen) {
      this.name = name;
      entries.add(new MetaAnalysisEntry(stud, phen));
    }

    @Override
    public int hashCode() {
      return Objects.hash(name);
    }

    @Override
    public boolean equals(Object obj) {
      if (this == obj) {
        return true;
      }
      if (!(obj instanceof MetaAnalysis)) {
        return false;
      }
      MetaAnalysis other = (MetaAnalysis) obj;
      return Objects.equals(name, other.name);
    }

  }

  private void loadMetaAnalysisFile() {
    if (!Files.exists(metaDir + META_ANALYSIS_DEFINITIONS_FILE)) {
      if (studies.size() > 1) {
        log.report("No meta analysis file (meta.analysis) found; no meta analysis will be run.");
      }
      return;
    }
    String[] metaAnalysisDefs = HashVec.loadFileToStringArray(metaDir
                                                              + META_ANALYSIS_DEFINITIONS_FILE,
                                                              false, null, false);

    for (String def : metaAnalysisDefs) {
      MetaAnalysis ma = new MetaAnalysis();
      String[] studDef = def.split("\t");
      if (studDef.length == 0) {
        continue;
      } else if (studDef.length == 1) {
        // TODO log - cannot meta-analyze with only one study!
        continue;
      }
      String name = studDef[0];
      for (int i = 1; i < studDef.length; i++) {
        String s = studDef[i];
        String[] studPhen = s.split("\\|");
        if (studPhen.length != 2) {
          // TODO report malformed study entry
          continue;
        }
        String stud = studPhen[0];
        String phen = studPhen[1];

        List<Study> match = studies.stream().filter(st -> st.studyName.equalsIgnoreCase(stud))
                                   .collect(Collectors.toList());
        if (match.size() > 1) {
          // TODO log malformed
          continue;
        } else if (match.size() == 0) {
          // TODO log not found study name
          continue;
        }

        Study matchStudy = match.get(0);

        List<String> matchPheno = matchStudy.phenoFiles.stream().filter(pde -> {
          return pde.equalsIgnoreCase(phen) || pde.equalsIgnoreCase(phen + ".pheno");
        }).collect(Collectors.toList());

        if (matchPheno.size() > 1) {
          // TODO log malformed
          continue;
        } else if (matchPheno.size() == 0) {
          // TODO log not found study name
          continue;
        }

        ma.addEntry(name, matchStudy, matchPheno.get(0));
      }
      if (ma.entries.size() <= 1) {
        // TODO log malformed
        continue;
      }

      metaAnalyses.add(ma);
    }

    for (MetaFile mf : metaFiles) {
      for (Constraint constr : analysisConstraints) {
        summaryResults.put(constr, mf, new HashMap<>());
      }
    }

    // TODO log entries found
  }

  private void runMetaHitWindowsAndLoadData() {
    for (MetaFile mf : metaFiles) {
      String metaFile = mf.metaFile;
      String fullPath = metaDir + metaFile;

      // we expect no duplicate marker names in meta file!
      Map<String, Alleles> metaHitMkrAlleles = new HashMap<>();
      Map<String, int[]> metaHitMkrLocations = new HashMap<>();
      Map<String, String[]> mkrAliasLookup = new HashMap<>();
      for (Constraint constr : analysisConstraints) {
        FileColumn<String> mkrColumn = StandardFileColumns.snp("mkr");
        FileColumn<Byte> chrColumn = StandardFileColumns.chr("chr");
        FileColumn<Integer> posColumn = StandardFileColumns.pos("pos");
        FileColumn<String> a1Column = StandardFileColumns.a1("a1");
        FileColumn<String> a2Column = StandardFileColumns.a2("a2");

        try {
          Map<ImmutableList<Object>, DataLine> data = FileParserFactory.setup(fullPath, mkrColumn,
                                                                              chrColumn, posColumn,
                                                                              a1Column, a2Column)
                                                                       .build()
                                                                       .load(false, mkrColumn,
                                                                             a1Column, a2Column);
          data.entrySet().stream().forEach(e -> {
            String mkr = e.getValue().get(mkrColumn, null);
            String a1 = e.getValue().get(a1Column, null);
            String a2 = e.getValue().get(a2Column, null);
            if (mkr == null || a1 == null || a2 == null) return;
            int[] chrPos = new int[] {e.getValue().get(chrColumn, (byte) -1),
                                      e.getValue().get(posColumn, -1)};
            metaHitMkrAlleles.put(mkr, new Alleles(a1, a2));
            metaHitMkrLocations.put(mkr, chrPos);
            String alias1 = "chr" + e.getValue().get(chrColumn, (byte) -1) + ":"
                            + e.getValue().get(posColumn, -1)
                            + (e.getValue().has(a1Column)
                                                          ? (":" + e.getValue().get(a1Column, "")
                                                             + ":" + e.getValue().get(a2Column, ""))
                                                          : "");
            String alias2 = "chr" + e.getValue().get(chrColumn, (byte) -1) + ":"
                            + e.getValue().get(posColumn, -1)
                            + (e.getValue().has(a1Column)
                                                          ? (":" + e.getValue().get(a2Column, "")
                                                             + ":" + e.getValue().get(a1Column, ""))
                                                          : "");
            mkrAliasLookup.put(alias1, new String[] {mkr, a1, a2});
            mkrAliasLookup.put(alias2, new String[] {mkr, a1, a2});
          });
          mf.metaMarkerAliasLookup = mkrAliasLookup;
        } catch (IOException e) {
          log.reportIOException(fullPath);
        }

        if (metaHitMkrLocations.isEmpty()) {
          log.reportError(".meta file was empty for " + metaFile);
        } else {
          log.report("Using all " + metaHitMkrLocations.size() + " SNPs in .meta file "
                     + mf.metaRoot);
        }

        int metaCount = Files.countLines(fullPath, 0);

        if (metaCount > 1000 && runMetaHW) {
          log.report("Found more than " + 1000 + " SNPs in meta file " + metaFile
                     + "; running HitWindows to trim...");
          metaHitMkrLocations.clear();
          String[][] results = HitWindows.determine(fullPath, constr.indexThreshold,
                                                    constr.windowMinSizePerSide,
                                                    constr.windowExtensionThreshold,
                                                    DEFAULT_ADDL_ANNOT_VAR_NAMES, new Logger());
          if (results == null) {
            log.reportError("HitWindows result was null for " + metaFile + ". Using all SNPs");
          } else {
            log.report("Found " + results.length + " hit windows");
            for (int i = 1; i < results.length; i++) {
              String[] result = results[i];
              try {
                metaHitMkrLocations.put(result[1], new int[] {Integer.parseInt(result[2]),
                                                              Integer.parseInt(result[3])});
              } catch (NumberFormatException nfe) {
                log.reportError("Failed to parse position for result " + ArrayUtils.toStr(result));
                metaHitMkrLocations.put(result[1], new int[] {-1, -1});
              }
            }
          }
        }

        // uncomment to use all markers in dataFile

        if (!metaHitMkrLocations.isEmpty()) {
          // read betas and freqs for hitwindow markers
          HashMap<String, MetaMarker> metaDataMarkers = new HashMap<>();
          Set<String> mkrs = new HashSet<>(metaHitMkrLocations.keySet());
          for (String hitMkr : mkrs) {
            if (mf.metaMarkers.containsKey(hitMkr)) {
              MetaMarker mm = mf.metaMarkers.get(hitMkr);
              if (!Double.isNaN(mm.beta) && !Double.isNaN(mm.freq) && !Double.isNaN(mm.pval)) {
                metaDataMarkers.put(hitMkr, mm);
              } else {
                metaHitMkrLocations.remove(hitMkr);
              }
            } else {
              metaHitMkrLocations.remove(hitMkr);
            }
          }

          double dataScore1 = getBetaFreqScore(metaDataMarkers);
          double dataScore2 = getChiDistRevScore(metaDataMarkers);

          // cross-ref PLINK markers
          for (Study study : studies) {
            // data sources MAY have duplicate markers!
            study.loadDataSources(mf, constr, metaHitMkrLocations);

            HashMap<String, MetaMarker> metaMarkerInDataSubset = new HashMap<>();
            Set<String> markerNames = new HashSet<>();
            markerNames.addAll(metaHitMkrLocations.keySet());
            markerNames.addAll(study.data.values().stream().map(d -> d.getMarkerSet())
                                         .flatMap(sms -> Arrays.stream(sms.getMarkerNames()))
                                         .collect(Collectors.toSet()));
            Multimap<String, String[]> bimMkrSet = study.retrieveMarkers(mf, constr, markerNames);

            // pull betas and freqs for union markers
            bimMkrSet.entries().stream().filter(Predicates.or(e -> {
              // filter on marker names
              return metaHitMkrLocations.containsKey(e.getKey())
                     || mf.metaMarkerAliasLookup.containsKey(e.getKey());
            }, e -> {
              // filter on marker alleles
              String[] dataAlleles = e.getValue();

              Alleles metaAlleles = metaHitMkrAlleles.get(e.getKey());
              if (metaAlleles == null && mf.metaMarkerAliasLookup.containsKey(e.getKey())) {
                metaAlleles = metaHitMkrAlleles.get(mf.metaMarkerAliasLookup.get(e.getKey())[0]);
              }
              if (metaAlleles == null) return false;

              return (dataAlleles[0].equalsIgnoreCase(metaAlleles.a1)
                      && dataAlleles[1].equalsIgnoreCase(metaAlleles.a2))
                     || (dataAlleles[0].equalsIgnoreCase(metaAlleles.a2)
                         && dataAlleles[1].equalsIgnoreCase(metaAlleles.a1));
            })).forEach(e -> {
              String dataMkr = e.getKey();

              if (!metaHitMkrLocations.containsKey(dataMkr)) {
                dataMkr = mf.metaMarkerAliasLookup.get(dataMkr)[0];
              }

              // because of the filter() above, we know metaAlleles is valid here
              metaMarkerInDataSubset.put(dataMkr, metaDataMarkers.get(dataMkr));
            });

            // apply equation and set value for overall results
            double bimScore1 = getBetaFreqScore(metaMarkerInDataSubset);
            double bimScore2 = getChiDistRevScore(metaMarkerInDataSubset);

            study.scores.put(constr, mf,
                             new double[] {bimScore1 / dataScore1, bimScore2 / dataScore2});
          }
        }
      }
    }
  }

  private void createAffectedPhenoFiles(Study study) {
    String affFile = study.studyDir + "AFFECTED.pheno";

    // skip if we don't have PHENO data
    if (!study.isPlinkData()) {
      return;
    }
    // if affected.pheno file already exists, skip
    if (study.phenoFiles.contains("AFFECTED.pheno")) {
      return;
    }
    // if any of the data folders have been created, skip creating affected.pheno file
    for (MetaFile mf : metaFiles) {
      if ((new File(study.studyDir + mf.metaRoot + "/")).exists()) {
        return;
      }
    }

    BufferedReader reader;
    // 0 - fid
    // 1 - iid
    // 2
    // 3
    // 4 - sex
    // 5 - pheno4

    ArrayList<String> fam = new ArrayList<>();
    Map<String, String> sex = new HashMap<>();
    Map<String, String> pheno = new HashMap<>();
    String temp;
    boolean override = false;
    try {
      // possible for `isPlinkData()` to return true if there are multiple PLINK data sets
      for (int i = 0; i < study.dataSources.size(); i++) {
        reader = Files.getAppropriateReader(study.dataSources.get(i).idFile);
        while ((temp = reader.readLine()) != null) {
          String[] line = temp.split(PSF.Regex.GREEDY_WHITESPACE);
          if (!ext.isMissingValue(line[5]) && ext.isValidDouble(line[5])) {
            String pheS = line[5];
            String sexS = "" + (ext.isMissingValue(line[4]) ? "."
                                                            : -1 * (Integer.parseInt(line[4]) - 2));
            fam.add(line[0] + "\t" + line[1]);
            sex.put(line[0] + "\t" + line[1], sexS);
            String prev = pheno.put(line[0] + "\t" + line[1], pheS);
            if (prev != null && !ext.isMissingValue(prev)) {
              log.reportError("Duplicate ID found in Plink data sets for study '" + study.studyName
                              + "'.  Please reconcile or combine data sets externally and try again.  GeneScorePipeline will continue without Plink phenotype analysis.");
              override = true;
              break;
            }
          }
        }
        if (override) break;
      }
    } catch (NumberFormatException e) {
      log.reportException(e);
      e.printStackTrace();
    } catch (IOException e) {
      log.reportException(e);
      e.printStackTrace();
    }

    String[] unique = ArrayUtils.unique(pheno.values().toArray(new String[] {}));
    if (unique.length == 1) {
      log.reportWarning("Error - no variance in pheno data from .fam file(s) for study '"
                        + study.studyName + "'.  Plink phenotype analysis will be skipped.");
      override = true;
    }
    if (override) {
      // either duplicate pheno entry or no variance, either way skip creating file
      return;
    }
    String header = "FID\tIID\tPHENO\tMALE";
    PrintWriter writer = Files.getAppropriateWriter(affFile);
    writer.println(header);
    for (String line : fam) {
      writer.println(line + "\t" + pheno.get(line) + "\t" + sex.get(line));
    }
    writer.flush();
    writer.close();

    study.phenoFiles.add("AFFECTED.pheno");
  }

  private void loadPhenoFiles(Study study) {
    for (String pheno : study.phenoFiles) {
      PhenoData pd = new PhenoData();
      pd.phenoName = pheno;

      try {
        BufferedReader reader = Files.getAppropriateReader(study.studyDir + pheno);
        String[] header = reader.readLine().trim().split("\t", -1);
        // fid == header[0]
        // iid == header[1]
        // depvar = header[2];
        if (header.length < 3) {
          throw new IllegalArgumentException("Couldn't parse three tab-delimited columns from pheno file: "
                                             + pheno);
        }
        ArrayList<String> covars = new ArrayList<>();
        for (int i = 3; i < header.length; i++) {
          covars.add(header[i]);
        }
        pd.covars.addAll(covars);
        String temp = reader.readLine();
        long count = 0;
        do {
          temp = temp.trim();
          count++;
          if (temp.isEmpty()) continue;

          String[] line = temp.split("\t", -1);

          try {
            if (!ext.isMissingValue(line[2])) {
              String fid = line[0];
              String iid = line[1];
              double depvar = Double.parseDouble(line[2]);
              PhenoIndiv pi = new PhenoIndiv(fid, iid, depvar);
              for (int i = 3; i < line.length; i++) {
                if (!ext.isMissingValue(line[i])) {
                  pi.addCovar(header[i], Double.parseDouble(line[i]));
                }
              }
              pd.indivs.put(pi.getFid() + "\t" + pi.getIid(), pi);
            }
          } catch (ArrayIndexOutOfBoundsException e) {
            throw new IllegalArgumentException("Error - found " + line.length
                                               + " columns (expected " + header.length
                                               + ") on line " + count, e);
          }
        } while ((temp = reader.readLine()) != null);

        reader.close();
      } catch (IOException e) {
        e.printStackTrace();
      }

      study.phenoData.put(pheno, pd);
    }
  }

  private double getBetaFreqScore(HashMap<String, MetaMarker> markerMap) {
    double sum = 0.0;
    // (Beta^2 * 2 * MAF (1-MAF))
    for (java.util.Map.Entry<String, MetaMarker> entry : markerMap.entrySet()) {
      MetaMarker mm = entry.getValue();
      double score = 2 * (mm.beta * mm.beta) * mm.freq * (1 - mm.freq);
      sum += score;
    }
    return sum;
  }

  private double getChiDistRevScore(HashMap<String, MetaMarker> markerMap) {
    double sum = 0.0;
    for (java.util.Map.Entry<String, MetaMarker> entry : markerMap.entrySet()) {
      sum += ProbDist.ChiDistReverse(entry.getValue().pval, 1);
    }
    return sum;
  }

  public void runPipeline() {
    log.report("Processing study data [" + studies.size() + " total]:");
    for (Study study : studies) {
      createAffectedPhenoFiles(study);
      loadPhenoFiles(study);
      processStudy(study);
    }
    processMetaAnalyses();
    writeResults();
    writeAndRunMR();
    log.report("Processing Complete!");
  }

  private void processStudy(Study study) {
    try {
      createFolders(study);
      crossFilterMarkerData(study);
      runHitWindows(study);
      extractHitMarkerData(study);
      runScore(study);
      runRegression(study);
      study.dumpStudy();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  private void createFolders(Study study) {
    for (MetaFile mf : metaFiles) {
      String dataFolder = study.studyDir + mf.metaRoot + "/";
      for (Constraint constr : analysisConstraints) {
        String constraintFolder = dataFolder + constr.analysisString + "/";
        File f = new File(constraintFolder);
        if (!(f.exists())) {
          f.mkdirs();
        }
      }
    }
  }

  private void crossFilterMarkerData(Study study) throws IOException {
    for (MetaFile mf : metaFiles) {
      String metaFile = mf.metaFile;
      String dataFile = mf.metaRoot;
      for (Constraint constraint : analysisConstraints) {
        String key = constraint.analysisString;
        String crossFilterFile = study.studyDir + dataFile + "/" + key + "/"
                                 + CROSS_FILTERED_DATAFILE;
        final boolean exists = new File(crossFilterFile).exists();
        if (exists && overwrite) {
          log.report("Cross-filtered data file already exists, but overwrite flag was set; recreating ...");
        }
        if (exists && !overwrite) {
          log.report("Cross-filtered data file already exists! [ --> '" + crossFilterFile + "']");
          study.hitSnpCounts.put(constraint, mf, Files.countLines(crossFilterFile, 1));
          continue;
        }
        if (study.data.get(constraint, mf).isEmpty()) {
          study.hitSnpCounts.put(constraint, mf, 0);
          continue;
        }

        log.report("Cross-filtering data and .meta files [ --> '" + crossFilterFile + "']");
        HashMap<String, GenomicPosition> mkrsMeta = new HashMap<>();
        Map<String, String[]> mkrAlleles = new HashMap<>();
        SnpMarkerSet markerSet = study.data.get(constraint, mf).getMarkerSet();
        List<Integer> inval = new ArrayList<>();
        List<Integer> ambig = new ArrayList<>();

        String[] mkrNames = markerSet.getMarkerNames();
        String[][] alleles = markerSet.getAlleles();
        int[][] chrPos = markerSet.getChrAndPositionsAsInts();

        for (int i = 0; i < mkrNames.length; i++) {
          String a1 = (alleles[i][0]).toUpperCase();
          String a2 = (alleles[i][1]).toUpperCase();
          boolean validAlleles = (Sequence.validAllele(a1) && Sequence.validAllele(a2));
          if (validAlleles) {
            mkrsMeta.put(mf.metaMarkerAliasLookup.containsKey(mkrNames[i]) ? mf.metaMarkerAliasLookup.get(mkrNames[i])[0]
                                                                           : mkrNames[i],
                         new GenomicPosition((byte) chrPos[i][0], chrPos[i][1]));
            mkrAlleles.put(mf.metaMarkerAliasLookup.containsKey(mkrNames[i]) ? mf.metaMarkerAliasLookup.get(mkrNames[i])[0]
                                                                             : mkrNames[i],
                           new String[] {a1, a2});
            if (a1.equals(Sequence.flip(a2))) {
              ambig.add(i);
            }
          } else {
            inval.add(i);
          }
        }

        if (ambig.size() > 0) {
          log.report("Found " + ambig.size() + " ambiguous markers out of " + mkrNames.length
                     + " (either A/T or G/C allele pairs); make sure you have everything on the same strand, or you many run into problems!");
          if (ambig.size() < 10) {
            StringBuilder builder = new StringBuilder(mkrNames[ambig.get(0)]);
            for (int i = 1; i < ambig.size(); i++) {
              builder.append(", ").append(mkrNames[ambig.get(i)]);
            }
            log.report("Ambiguous markers: " + builder.toString());
          } else {
            PrintWriter writer = Files.getAppropriateWriter(metaDir + "ambiguous.txt");
            for (Integer i : ambig) {
              writer.println(mkrNames[i] + "\t" + chrPos[i][0] + "\t" + chrPos[i][1] + "\t"
                             + alleles[i][0] + "\t" + alleles[i][1]);
            }
            writer.close();
          }
        }
        if (inval.size() > 0) {
          log.report("Found " + inval.size() + " invalid markers out of " + mkrNames.length + ".");
          if (inval.size() < 10) {
            StringBuilder builder = new StringBuilder(inval.get(0));
            for (int i = 1; i < inval.size(); i++) {
              builder.append(", ").append(inval.get(i));
            }
            log.report("Invalid markers: " + builder.toString());
          } else {
            PrintWriter writer = Files.getAppropriateWriter(metaDir + "invalid.txt");
            for (Integer i : inval) {
              writer.println(mkrNames[i] + "\t" + chrPos[i][0] + "\t" + chrPos[i][1] + "\t"
                             + alleles[i][0] + "\t" + alleles[i][1]);
            }
            writer.close();
          }
        }

        final String outDelim = "\t";

        final AliasedFileColumn markerCol = StandardFileColumns.snp("MarkerName");
        final FileColumn<Byte> chrLinkedColumn = new AbstractFileColumn<Byte>("Chr") {

          @Override
          public Byte getValue(String[] line) throws ParseFailureException {
            String mkr = markerCol.getValue(line);
            if (mkr == null) {
              throw new ParseFailureException("Marker not found in line: " + line.toString());
            }
            GenomicPosition gp = mkrsMeta.get(mkr);
            if (gp == null) {
              throw new ParseFailureException("Marker not found in data: " + mkr);
            }
            return gp.getChr();
          }
        };
        final FileColumn<Integer> posLinkedColumn = new AbstractFileColumn<Integer>("Position") {

          @Override
          public Integer getValue(String[] line) throws ParseFailureException {
            String mkr = markerCol.getValue(line);
            if (mkr == null) {
              throw new ParseFailureException("Marker not found in line: " + line.toString());
            }
            GenomicPosition gp = mkrsMeta.get(mkr);
            if (gp == null) {
              throw new ParseFailureException("Marker not found in data: " + mkr);
            }
            return mkrsMeta.get(mkr).getPosition();
          }
        };

        final AliasedFileColumn a1Column = StandardFileColumns.a1("a1");
        final AliasedFileColumn a2Column = StandardFileColumns.a2("a2");
        final IndexedFileColumn<Double> pColumn = StandardFileColumns.pVal("p");
        final FileColumn<String> remainingColumns = StandardFileColumns.allExcept(outDelim,
                                                                                  markerCol,
                                                                                  pColumn,
                                                                                  StandardFileColumns.chr("excludedChr"),
                                                                                  StandardFileColumns.pos("excludedPos"));
        final ColumnFilter mkrsBimFilter = new AbstractColumnFilter(markerCol) {

          @Override
          public boolean filter(DataLine values) {
            String mkr = values.get(markerCol, null);
            return mkrsMeta.containsKey(mkr)
                   || (mf.metaMarkerAliasLookup.containsKey(mkr)
                       && mkrsMeta.containsKey(mf.metaMarkerAliasLookup.get(mkr)[0]));
          }
        };
        final ColumnFilter mkrsAlleleFilter = new AbstractColumnFilter(markerCol, a1Column,
                                                                       a2Column) {

          @Override
          public boolean filter(DataLine values) {
            String mkr = values.get(markerCol, null);
            String[] dataAlleles = mkrAlleles.get(mkr);
            String[] metaAlleles = new String[] {values.get(a1Column, "NA"),
                                                 values.get(a2Column, "NA")};
            CONFIG config = StrandOps.determineStrandConfig(dataAlleles, metaAlleles);
            AlleleOrder alleleOrder = config.getAlleleOrder();
            if (config == CONFIG.AMBIGUOUS) {
              boolean a1 = dataAlleles[0].equalsIgnoreCase(metaAlleles[0])
                           || dataAlleles[0].equalsIgnoreCase(metaAlleles[1]);
              boolean a2 = dataAlleles[1].equalsIgnoreCase(metaAlleles[0])
                           || dataAlleles[1].equalsIgnoreCase(metaAlleles[1]);
              return a1 && a2;
            } else if (alleleOrder == AlleleOrder.SAME || alleleOrder == AlleleOrder.OPPOSITE)
              return true;
            else {
              Joiner alleleJoiner = Joiner.on('/');
              log.reportError("Alleles in study (" + alleleJoiner.join(dataAlleles)
                              + ") do not match source alleles (" + alleleJoiner.join(metaAlleles)
                              + ") for " + mkr);
              return false;
            }
          }
        };
        final ColumnFilter pValThreshold = new DoubleFilter(pColumn, COMPARISON.LT,
                                                            constraint.indexThreshold);
        final ColumnFilter pValMissing = new AbstractColumnFilter(pColumn) {

          @Override
          public boolean filter(DataLine values) {
            return !values.has(pColumn);
          }
        };
        final ColumnFilter pValFilter = ColumnFilters.or(pValMissing, pValThreshold);
        try (FileParser crossFilterParser = FileParserFactory.setup(metaDir + metaFile, markerCol,
                                                                    chrLinkedColumn,
                                                                    posLinkedColumn,
                                                                    remainingColumns)
                                                             .optionalColumns(pColumn)
                                                             .filter(mkrsBimFilter, false, true)
                                                             .filter(mkrsAlleleFilter, false, true)
                                                             .filter(pValFilter, false, true)
                                                             .build()) {

          int hitCount = crossFilterParser.parseToFile(crossFilterFile, outDelim);
          int count;

          count = crossFilterParser.getFilteredCount(pValFilter);
          log.report("Dropped " + count + " snps for p-value > " + constraint.indexThreshold
                     + (count > 0 ? (": {" + crossFilterParser.getFilteredData(pValFilter).stream()
                                                              .map(d -> d.getUnsafe(markerCol))
                                                              .collect(Collectors.joining(", "))
                                     + "}.")
                                  : "."));

          count = crossFilterParser.getFilteredCount(mkrsAlleleFilter);
          log.report("Dropped " + count + " snps for mismatched alleles"
                     + (count > 0 ? (": {"
                                     + crossFilterParser.getFilteredData(mkrsAlleleFilter).stream()
                                                        .map(d -> d.getUnsafe(markerCol))
                                                        .collect(Collectors.joining(", "))
                                     + "}.")
                                  : "."));

          count = crossFilterParser.getFilteredCount(mkrsBimFilter);
          log.report("Dropped " + count + " snps for lacking data"
                     + (count > 0 ? (": {"
                                     + crossFilterParser.getFilteredData(mkrsBimFilter).stream()
                                                        .map(d -> d.getUnsafe(markerCol))
                                                        .collect(Collectors.joining(", "))
                                     + "}.")
                                  : "."));

          study.hitSnpCounts.put(constraint, mf, hitCount);
        }
      }
    }
  }

  private void runHitWindows(Study study) {
    for (MetaFile mf : metaFiles) {

      for (Constraint constr : analysisConstraints) {
        String analysisKey = constr.analysisString;
        File prefDir = new File(getDirPath(study, mf, constr));
        String crossFilterFile = prefDir + "/" + CROSS_FILTERED_DATAFILE;
        String hitsFile = prefDir + FILE_PREFIX_HITS + analysisKey + FILE_EXT_OUT;
        final boolean exists = (new File(hitsFile)).exists();
        if (exists && overwrite) {
          log.report("Hit window analysis file already exists, but overwrite flag was set; recreating ... ");
        }
        if (exists && !overwrite) {
          log.report("Hit window analysis file already exists! [ --> '" + hitsFile + "']");
          study.hitWindowCnts.put(constr, mf, Files.countLines(hitsFile, 1));
          continue;
        }
        if (study.data.get(constr, mf).isEmpty()) {
          continue;
        }
        if (!Files.exists(crossFilterFile, true)) {
          log.reportError("cross filter data file (" + crossFilterFile
                          + ") is either missing or empty; cannot run HitWindows on an empty file!");
          continue;
        }

        log.report("Running hit window analysis [ --> '" + hitsFile + "']");
        String[][] results = HitWindows.determine(crossFilterFile, constr.indexThreshold,
                                                  constr.windowMinSizePerSide,
                                                  constr.windowExtensionThreshold,
                                                  DEFAULT_ADDL_ANNOT_VAR_NAMES, log);
        if (results == null) {
          log.reportError("HitWindows result from " + crossFilterFile + " was null");
        } else {
          int hits = results.length - 1; // Don't count header
          log.report("Found " + hits + " hit windows");
          Files.writeMatrix(results, hitsFile, "\t");
          study.hitWindowCnts.put(constr, mf, hits);
        }
      }
    }
  }

  private void extractHitMarkerData(Study study) {
    for (MetaFile mf : metaFiles) {

      for (Constraint constr : analysisConstraints) {
        if (study.data.get(constr, mf).isEmpty()) {
          study.markerData.put(constr, mf, new HashMap<>());
          continue;
        }

        File prefDir = new File(getDirPath(study, mf, constr));

        String crossFilterFile = prefDir + "/" + CROSS_FILTERED_DATAFILE;
        String hitsFile = prefDir + FILE_PREFIX_HITS + constr.analysisString + FILE_EXT_OUT;
        String[] hitMarkers = HashVec.loadFileToStringArray(hitsFile, true,
                                                            new int[] {hitsMkrIndex}, false);
        HashSet<String> hitMrkSet = new HashSet<>();
        for (String mkr : hitMarkers) {
          hitMrkSet.add(mkr);
        }

        final FileColumn<String> markerCol = StandardFileColumns.snp(MARKER_COL_NAME);
        final FileColumn<?> effectAlleleCol = StandardFileColumns.a1(EFFECT_ALLELE_COL_NAME);
        final FileColumn<?> nonEffectAlleleCol = StandardFileColumns.a2(NON_EFFECT_ALLELE_COL_NAME);
        final FileColumn<?> betaCol = StandardFileColumns.beta(BETA_COL_NAME);
        ColumnFilter hitMarkerFilter = new AbstractColumnFilter(markerCol) {

          @Override
          public boolean filter(DataLine values) {
            return hitMrkSet.contains(values.get(markerCol, null));
          }
        };
        try (FileParser crossFilterParser = FileParserFactory.setup(crossFilterFile, markerCol)
                                                             .optionalColumns(effectAlleleCol,
                                                                              nonEffectAlleleCol,
                                                                              betaCol)
                                                             .filter(hitMarkerFilter).build()) {
          String subsetFile = prefDir + "/subsetData_" + constr.analysisString + ".xln";
          Map<String, HitMarker> dataList = crossFilterParser.parseToFileAndLoad(subsetFile, "\t",
                                                                                 true, markerCol)
                                                             .entrySet().stream()
                                                             .collect(Collectors.toMap(Map.Entry::getKey,
                                                                                       e -> formHitMarker(e.getValue())));
          study.markerData.put(constr, mf, dataList);
        } catch (IOException e) {
          log.reportIOException(crossFilterFile);
        }
      }
    }
  }

  private static final HitMarker formHitMarker(DataLine hitMarkerDataLine) {
    return new HitMarker(hitMarkerDataLine.get(StandardFileColumns.a1(EFFECT_ALLELE_COL_NAME),
                                               "NA"),
                         hitMarkerDataLine.get(StandardFileColumns.a2(NON_EFFECT_ALLELE_COL_NAME),
                                               "NA"),
                         hitMarkerDataLine.get(StandardFileColumns.beta(BETA_COL_NAME), null));
  }

  private static final String MARKER_REGRESSION_PREFIX = "marker_results_";
  private static final String MARKER_REGRESSION_EXTEN = FILE_EXT_OUT;
  private static final String MENDELIAN_RANDOMIZATION_PREFIX = "MendelianRandomization_";
  private static final String TWO_SAMPLE_PREFIX = "TwoSampleMR_";
  private static final String TWO_SAMPLE_HDR = "SNP\teffect_allele\tother_allele\teaf\tbeta\tse\tP-value\tPhenotype";
  private static final String MENDELIAN_RANDOMIZATION_EXTEN = ".csv";
  private static final String REGRESSION_BETA_FILE = "markers_beta.out";

  private void runScore(Study study) {
    for (MetaFile mf : metaFiles) {
      String dataFile = mf.metaRoot;
      for (Constraint constr : analysisConstraints) {
        String analysisKey = constr.analysisString;

        log.report("Setting up for PRS computation for {" + dataFile + "\t" + analysisKey + "}");

        File prefDir = new File(getDirPath(study, mf, constr));
        if (!prefDir.exists()) {
          log.report("Error - no subfolder for '" + analysisKey + "' analysis");
          continue;
        }
        DosageData data = study.data.get(constr, mf);
        if (data.isEmpty()) {
          log.reportError("Dosage data was empty for {" + dataFile + "\t" + analysisKey + "}");
          continue;
        }

        // this data has been filtered already to remove duplicate markers
        // therefore it is safe to key on marker name
        // however we still have to match up against data alleles
        Map<String, HitMarker> hitMarkerData = study.markerData.get(constr, mf);

        String[][] ids = data.getIds();
        float[][] dose = data.getDosageValues();
        if (dose == null) {
          if (data.hasGenotypeProbabilities()) {
            data.computeDosageValues(log);
            dose = data.getDosageValues();
          }
        }

        if (dose == null || dose.length == 0) {
          log.reportError("No dosage data available for {" + dataFile + "\t" + analysisKey + "}");
          continue;
        }

        String[] markers = data.getMarkerSet().getMarkerNames();
        String[][] alleles = data.getMarkerSet().getAlleles();

        study.markerDosageIndices.put(constr, mf, HashBasedTable.create());
        study.matchedMarkerIndices.put(constr, mf, new HashMap<>());
        Table<String, Alleles, Integer> markerDosageIndices = study.markerDosageIndices.get(constr,
                                                                                            mf);
        Map<String, Integer> matchedMarkerIndices = study.matchedMarkerIndices.get(constr, mf);

        Map<String, Float> matchedMarkerFreqs = new HashMap<>();
        Map<String, Alleles> matchedMarkerAlleles = new HashMap<>();
        Map<String, AlleleOrder> matchedMarkerAlleleOrders = new HashMap<>();

        int index = 0;
        for (int m = 0; m < markers.length; m++) {
          String mkr = markers[m];
          String[] dataAlleles = alleles[m];
          Alleles dataAlleleObj = new Alleles(dataAlleles[0], dataAlleles[1]);
          markerDosageIndices.put(mkr, dataAlleleObj, m);
          HitMarker hitMarker = hitMarkerData.get(mkr);
          if (hitMarker == null) {
            if (mf.metaMarkerAliasLookup.containsKey(mkr)) {
              hitMarker = hitMarkerData.get(mf.metaMarkerAliasLookup.get(mkr)[0]);
            }
            if (hitMarker == null) {
              log.report("No HitMarker data available for " + mkr);
              continue;
            }
          }
          boolean allelesMatch = (hitMarker.effectAllele.equalsIgnoreCase(dataAlleles[0])
                                  && hitMarker.nonEffectAllele.equalsIgnoreCase(dataAlleles[1]))
                                 || (hitMarker.effectAllele.equalsIgnoreCase(dataAlleles[1])
                                     && hitMarker.nonEffectAllele.equalsIgnoreCase(dataAlleles[0]));
          if (!allelesMatch) {
            continue;
          }
          CONFIG config = determineAlleleConfig(alleles[m], hitMarker);
          AlleleOrder alleleOrder = config.getAlleleOrder();
          if (config == CONFIG.AMBIGUOUS) {
            alleleOrder = hitMarker.getEffectAllele()
                                   .compareToIgnoreCase(alleles[m][0]) == 0 ? StrandOps.AlleleOrder.SAME
                                                                            : StrandOps.AlleleOrder.OPPOSITE;
          } else if (!(alleleOrder.equals(AlleleOrder.SAME)
                       || alleleOrder.equals(AlleleOrder.OPPOSITE))) {
            Joiner alleleJoiner = Joiner.on('/');
            log.reportError("Alleles in study (" + alleleJoiner.join(alleles[m])
                            + ") do not match source alleles ("
                            + alleleJoiner.join(hitMarker.getEffectAllele() == null ? "NULL"
                                                                                    : hitMarker.getEffectAllele(),
                                                hitMarker.getNonEffectAllele() == null ? "NULL"
                                                                                       : hitMarker.getNonEffectAllele())
                            + ") for " + mkr);
            continue;
          }
          matchedMarkerAlleleOrders.put(mkr, alleleOrder);
          matchedMarkerIndices.put(mkr, index++);
          matchedMarkerAlleles.put(mkr, dataAlleleObj);
          int cnt = 0;
          float tot = 0;
          for (int i = 0; i < ids.length; i++) {
            if (!Float.isNaN(dose[m][i])) {
              tot += dose[m][i];
              cnt++;
            }
          }
          matchedMarkerFreqs.put(mkr, tot / cnt);
        }

        log.report("Computing PRS for {" + dataFile + "\t" + analysisKey + "} with "
                   + matchedMarkerIndices.size() + " markers and " + ids.length + " samples...");

        double[][] markerScores = new double[ids.length][];
        double[][] markerDosages = new double[ids.length][];
        study.matchedMarkerScores.put(constr, mf, markerScores);
        study.matchedMarkerDosages.put(constr, mf, markerDosages);

        for (int i = 0; i < ids.length; i++) {
          double scoreSum = 0;
          double cnt2 = 0;
          int cnt = 0;

          int numMkrs = matchedMarkerIndices.size();

          markerScores[i] = new double[numMkrs];
          markerDosages[i] = new double[numMkrs];

          for (Map.Entry<String, Integer> mkrIndexEntry : matchedMarkerIndices.entrySet()) {
            String mkr = mkrIndexEntry.getKey();
            Alleles mkrAlleles = matchedMarkerAlleles.get(mkr);
            int mkrDosageIndex = markerDosageIndices.get(mkr, mkrAlleles);
            int mkrMatchedIndex = mkrIndexEntry.getValue();
            float dosage = dose[mkrDosageIndex][i];
            double mkrDose = Double.NaN;
            boolean isNaN = Float.isNaN(dosage);
            float mkrFrq = matchedMarkerFreqs.get(mkr);
            boolean nanFrq = Float.isNaN(mkrFrq);
            AlleleOrder alleleOrder = matchedMarkerAlleleOrders.get(mkr);

            HitMarker hitMarker = hitMarkerData.get(mkr);
            if (hitMarker == null) {
              if (mf.metaMarkerAliasLookup.containsKey(mkr)) {
                mkr = mf.metaMarkerAliasLookup.get(mkr)[0];
                hitMarker = hitMarkerData.get(mkr);
              }
              if (hitMarker == null) {
                log.report("No HitMarker data available for " + mkr
                           + ((mkrIndexEntry.getKey().equals(mkr)) ? ""
                                                                   : (" (" + mkrIndexEntry.getKey()
                                                                      + ")")));
                markerScores[i][mkrMatchedIndex] = Double.NaN;
                markerDosages[i][mkrMatchedIndex] = Double.NaN;
                continue;
              }
            }
            double beta = hitMarker.getEffect();
            double mkrScr = Double.NaN;

            if (alleleOrder.equals(StrandOps.AlleleOrder.OPPOSITE)) {
              cnt += isNaN ? 0 : 1;
              cnt2 += isNaN ? 0 : (mkrDose = dosage);
              mkrScr = (isNaN ? (nanFrq ? 0 : mkrFrq) : dosage) * beta;
            } else if (alleleOrder.equals(StrandOps.AlleleOrder.SAME)) {
              cnt += isNaN ? 0 : 1;
              cnt2 += isNaN ? 0 : (mkrDose = (2.0 - dosage));
              mkrScr = ((2.0 - (isNaN ? (nanFrq ? 0 : mkrFrq) : dosage)) * beta);
            } else {
              throw new IllegalStateException("Mismatched alleles were not caught when cross-filtering");
            }
            scoreSum += mkrScr;

            markerScores[i][mkrMatchedIndex] = mkrScr;
            markerDosages[i][mkrMatchedIndex] = mkrDose;

          }

          double mkrRatio = cnt / (double) numMkrs;
          if (mkrRatio > minMissThresh) {
            study.indivMetaScores.get(constr, mf).put(ids[i][0] + "\t" + ids[i][1], scoreSum);
            study.indivScoreData.put(constr, ids[i][0] + "\t" + ids[i][1],
                                     new double[] {mkrRatio, 2 * cnt, cnt2});
          } else {
            log.reportWarning("Missingness threshold not met for ID (" + ids[i][0] + " | "
                              + ids[i][1] + ")");
          }
        }

      }
    }
  }

  private static CONFIG determineAlleleConfig(String[] alleles, HitMarker hitMarker) {
    CONFIG config = StrandOps.determineStrandConfig(new String[] {hitMarker.getEffectAllele(),
                                                                  hitMarker.getNonEffectAllele()},
                                                    alleles);
    return config;
  }

  private void runRegression(Study study) {
    for (MetaFile mf : metaFiles) {
      for (Constraint constr : analysisConstraints) {
        if (study.data.get(constr, mf).isEmpty()) {
          continue;
        }

        DosageData data = study.data.get(constr, mf);
        String[][] ids = data.getIds();
        TrioScoreTest trioTests = new TrioScoreTest();

        int sibs = 0;

        for (Trio trio : data.getTrios()) {
          String caseID = ids[trio.getChildIndex()][PSF.Plink.FAM_FID_INDEX] + "\t"
                          + ids[trio.getChildIndex()][PSF.Plink.FAM_IID_INDEX];
          if (ids[trio.getChildIndex()].length > PSF.Plink.FAM_AFF_INDEX
              && !PSF.Plink.affIsCase(ids[trio.getChildIndex()][PSF.Plink.FAM_AFF_INDEX])) {
            sibs++;
            continue;
          }
          String fatherID = ids[trio.getFatherIndex()][PSF.Plink.FAM_FID_INDEX] + "\t"
                            + ids[trio.getFatherIndex()][PSF.Plink.FAM_IID_INDEX];
          String motherID = ids[trio.getMotherIndex()][PSF.Plink.FAM_FID_INDEX] + "\t"
                            + ids[trio.getMotherIndex()][PSF.Plink.FAM_IID_INDEX];
          if (study.indivMetaScores.get(constr, mf).containsKey(caseID)
              && study.indivMetaScores.get(constr, mf).containsKey(fatherID)
              && study.indivMetaScores.get(constr, mf).containsKey(motherID)) {
            double caseScore = study.indivMetaScores.get(constr, mf).get(caseID);
            double fatherScore = study.indivMetaScores.get(constr, mf).get(fatherID);
            double motherScore = study.indivMetaScores.get(constr, mf).get(motherID);
            trioTests.addScore(caseScore, fatherScore, motherScore);
          }
          TrioScoreTest.Results trioTestResults = trioTests.runTests();
          if (trioTestResults != null) study.trioScoreTests.put(constr, mf, trioTestResults);
        }

        if (sibs > 0) {
          log.reportWarning(sibs + " trios (of " + data.getTrios().size()
                            + " total trios) were excluded from paired score analyses because the children were not coded as cases in the fam data.");
        }

        Set<String> mkrs = study.matchedMarkerIndices.get(constr, mf).keySet();
        Set<String> intersectMkrs = Sets.intersection(mkrs, mf.metaMarkers.keySet());
        List<String> markersInOrder = new ArrayList<>(intersectMkrs);
        final Map<String, Integer> map = study.matchedMarkerIndices.get(constr, mf);
        markersInOrder.sort((m1, m2) -> {
          return map.getOrDefault(m1, -1).compareTo(map.getOrDefault(m2, -1));
        });

        PrintWriter writer = Files.getAppropriateWriter(metaDir + "missingPhenos.id");
        for (int i = 0; i < study.phenoFiles.size(); i++) {
          PhenoData pd = study.phenoData.get(study.phenoFiles.get(i));
          RegressionResult rrResult = actualRegression(study.indivMetaScores.get(constr, mf),
                                                       writer, pd);
          study.regressions.get(constr, mf).put(pd, rrResult);

          double[][] mkrDosages = study.matchedMarkerDosages.get(constr, mf);
          List<String> foundMetaMarkers = new ArrayList<>();
          for (String marker : markersInOrder) {
            if (mf.metaMarkers.containsKey(marker)) { // may have been dropped by meta HitWindows
              foundMetaMarkers.add(marker);
            }
          }

          String[][] dataIds = study.data.get(constr, mf).getIds();
          Map<String, double[]> indivAllMkrs = new HashMap<>();
          for (String[] ind : dataIds) {
            indivAllMkrs.put(ind[0] + "\t" + ind[1], new double[mkrDosages[0].length]);
          }

          List<String> mkrNames = new ArrayList<>();
          for (int indMkr = 0; indMkr < foundMetaMarkers.size(); indMkr++) {
            String marker = foundMetaMarkers.get(indMkr);
            mkrNames.add(marker);
            int mkrInd = study.matchedMarkerIndices.get(constr, mf).getOrDefault(marker, -1);
            if (mkrInd == -1) {
              // TODO probably log?
              continue;
            }

            double sum = 0;
            Map<String, Double> dosages = new HashMap<>();
            for (int s = 0; s < ids.length; s++) {
              final String key = ids[s][0] + "\t" + ids[s][1];
              dosages.put(key, mkrDosages[s][mkrInd]);
              indivAllMkrs.get(key)[indMkr] = mkrDosages[s][mkrInd];
              if (pd.indivs.containsKey(key)) {
                sum += mkrDosages[s][mkrInd];
              }
            }

            study.markerMAC.get(constr, mf).put(pd, marker, sum);

            if (study.mkrPassesMAC(constr, mf, pd, marker)) {
              RegressionResult rrResultPerMkr = actualRegression(dosages, null, pd);
              study.markerRegressions.get(constr, mf).put(pd, marker, rrResultPerMkr);
            } else {
              log.report("Marker " + marker + " did not pass CMAC filter for " + pd.phenoName);
            }
          }

          // run regression with all markers from indivAllMarkers
          Map<String, RegressionResult> rrResultMulti = actualRegressionMultivar(study, constr, mf,
                                                                                 indivAllMkrs,
                                                                                 mkrNames, null,
                                                                                 pd);
          rrResultMulti.forEach((s, rr1) -> {
            study.multivarRegressions.get(constr, mf).put(pd, s, rr1);
          });
        }
      }
    }
  }

  private Map<String, RegressionResult> actualRegressionMultivar(Study study, Constraint constr,
                                                                 MetaFile mf,
                                                                 Map<String, double[]> scoreData,
                                                                 List<String> scoreNames,
                                                                 PrintWriter writer, PhenoData pd) {
    ArrayList<Double> depData = new ArrayList<>();
    ArrayList<double[]> baselineIndeps = new ArrayList<>();
    ArrayList<double[]> indepData = new ArrayList<>();
    MultiSet<String> invalids = new HashMultiSet<>();
    for (java.util.Map.Entry<String, PhenoIndiv> indiv : pd.indivs.entrySet()) {
      if (scoreData.containsKey(indiv.getKey())) {
        PhenoIndiv pdi = pd.indivs.get(indiv.getKey());
        double[] indivScoreData = scoreData.get(indiv.getKey());
        double[] baseData = new double[pd.covars.size() + covarData.size()];
        double[] covarData1 = new double[pd.covars.size() + indivScoreData.length
                                         + covarData.size()];
        System.arraycopy(indivScoreData, 0, covarData1, 0, indivScoreData.length);
        boolean validCovars = true;
        for (int k = 0; k < pd.covars.size(); k++) {
          Double d = pdi.getCovars().get(pd.covars.get(k));
          if (d == null) {
            if (writer != null) {
              writer.println(indiv.getKey() + "\t" + pd.covars.get(k));
            }
            validCovars = false;
            invalids.add(pd.covars.get(k));
          } else {
            baseData[k] = d.doubleValue();
            covarData1[indivScoreData.length + k] = d.doubleValue();
          }
        }
        for (int k = 0; k < covarOrder.size(); k++) {
          String covarKey = covarOrder.get(k);
          Double d = covarData.get(covarKey).get(indiv.getKey());
          if (d == null || d.isNaN()) {
            if (writer != null) {
              writer.println(indiv.getKey() + "\t" + covarKey);
            }
            validCovars = false;
            invalids.add(covarKey);
          } else {
            baseData[pd.covars.size() + k] = d.doubleValue();
            covarData1[indivScoreData.length + pd.covars.size() + k] = d.doubleValue();
          }
        }
        if (validCovars) {
          depData.add(pdi.getDepvar());
          baselineIndeps.add(baseData);
          indepData.add(covarData1);
        }
      }
    }
    if (invalids.uniqueSet().size() > 0) {
      for (String k : invalids.uniqueSet()) {
        log.reportWarning("Missing " + invalids.getCount(k) + " for covariate " + k + ".");
      }
    }

    Set<String> removedMkrs = new HashSet<>();
    List<Integer> removeCols = new ArrayList<>();
    for (int i = 0; i < scoreNames.size(); i++) {
      if (!study.mkrPassesMAC(constr, mf, pd, scoreNames.get(i))) {
        removedMkrs.add(scoreNames.get(i));
        removeCols.add(i);
      }
    }
    if (removeCols.size() > 0) {
      log.reportWarning("Removing " + removeCols.size() + " low-variance independent variables ["
                        + ArrayUtils.toStr(removedMkrs, ", ") + "] for multivariate regression for "
                        + pd.phenoName + " ...");
      removeCols.sort(Integer::compare);
      for (int i = removeCols.size() - 1; i >= 0; i--) {
        for (int r = 0; r < indepData.size(); r++) {
          double[] newCov = new double[indepData.get(r).length - 1];
          int indx = 0;
          for (int c = 0; c < indepData.get(r).length; c++) {
            if (c == removeCols.get(i)) continue;
            newCov[indx++] = indepData.get(r)[c];
          }
          indepData.set(r, newCov);
        }
      }
    }

    double[][] baseCovars = new double[baselineIndeps.size()][];
    double[][] covars = new double[indepData.size()][];
    for (int k = 0; k < covars.length; k++) {
      covars[k] = indepData.get(k);
      baseCovars[k] = baselineIndeps.get(k);
    }

    int cases = 0;
    int controls = 0;
    Multiset<Double> phenos = ImmutableMultiset.copyOf(depData);
    if (phenos.entrySet().size() == 2) {
      if (phenos.elementSet().equals(ImmutableSet.of(Double.valueOf(0.0), Double.valueOf(1.0)))) {
        cases = phenos.count(Double.valueOf(1.0));
        controls = phenos.count(Double.valueOf(0.0));
      } else if (phenos.elementSet()
                       .equals(ImmutableSet.of(Double.valueOf(1.0), Double.valueOf(2.0)))) {
        cases = phenos.count(Double.valueOf(2.0));
        controls = phenos.count(Double.valueOf(1.0));
      } else {
        log.reportError("Unrecognized case/control designations, cannot count cases and controls");
      }
    }
    RegressionModel baseModel = RegressionModel.determineAppropriate(Doubles.toArray(depData),
                                                                     baseCovars, false, true);
    RegressionModel model = RegressionModel.determineAppropriate(Doubles.toArray(depData), covars,
                                                                 false, true);
    List<String> varNameList = new ArrayList<>();
    for (String v : scoreNames) {
      varNameList.add(v);
    }
    for (String v : pd.covars) {
      varNameList.add(v);
    }
    for (String v : covarOrder) {
      varNameList.add(v);
    }
    for (int i = removeCols.size() - 1; i >= 0; i--) {
      varNameList.remove((int) removeCols.get(i));
    }

    String[] varNames = varNameList.toArray(new String[0]);
    model.setVarNames(varNames);

    Map<String, RegressionResult> results = new HashMap<>();

    log.report("Found " + model.getVarNames().length + " variables for regression ("
               + ArrayUtils.toStr(model.getVarNames(), ",") + ")");

    for (String var : scoreNames) {
      if (removedMkrs.contains(var)) {
        continue;
      }
      RegressionResult.Builder rr = RegressionResult.builder();
      rr.setNum(depData.size());
      rr.setnCases(cases);
      rr.setnControls(controls);
      if (model.analysisFailed()) {
        rr.setBaseRSq(baseModel.analysisFailed() ? Double.NaN : baseModel.getRsquare());
        rr.setBeta(Double.NaN);
        rr.setSe(Double.NaN);
        rr.setRsq(Double.NaN);
        rr.setPval(Double.NaN);
        rr.setLogistic(model.isLogistic());
      } else {
        int ind = -1;
        for (int l = 0; l < model.getVarNames().length; l++) {
          if (var.equals(model.getVarNames()[l])) {
            ind = l;
            break;
          }
        }
        if (ind == -1) {
          rr.setBaseRSq(baseModel.analysisFailed() ? Double.NaN : baseModel.getRsquare());
          rr.setBeta(Double.NaN);
          rr.setSe(Double.NaN);
          rr.setRsq(model.getRsquare());
          rr.setPval(Double.NaN);
          rr.setLogistic(model.isLogistic());
          rr.setStats(Double.NaN);
        } else {
          rr.setBaseRSq(baseModel.analysisFailed() ? Double.NaN : baseModel.getRsquare());
          rr.setBeta(model.getBetas()[ind]);
          rr.setSe(model.getSEofBs()[ind]);
          rr.setRsq(model.getRsquare());
          rr.setPval(model.getSigs()[ind]);
          rr.setLogistic(model.isLogistic());
          rr.setStats(model.getStats()[ind]);
        }
      }
      rr.setDroppedSnps(removeCols.size());
      RegressionResult rrResult = rr.build();
      results.put(var, rrResult);
    }
    return results;
  }

  private RegressionResult actualRegression(Map<String, Double> scoreData, PrintWriter writer,
                                            PhenoData pd) {
    ArrayList<Double> depData = new ArrayList<>();
    ArrayList<double[]> baselineIndeps = new ArrayList<>();
    ArrayList<double[]> indepData = new ArrayList<>();
    MultiSet<String> invalids = new HashMultiSet<>();
    for (java.util.Map.Entry<String, PhenoIndiv> indiv : pd.indivs.entrySet()) {
      if (scoreData != null && scoreData.containsKey(indiv.getKey())) {
        PhenoIndiv pdi = pd.indivs.get(indiv.getKey());
        double[] baseData = new double[pd.covars.size() + covarData.size()];
        double[] covarData1 = new double[pd.covars.size() + 1 + covarData.size()];
        covarData1[0] = scoreData.get(indiv.getKey());
        boolean validCovars = true;
        for (int k = 1; k < pd.covars.size() + 1; k++) {
          Double d = pdi.getCovars().get(pd.covars.get(k - 1));
          if (d == null) {
            if (writer != null) {
              writer.println(indiv.getKey() + "\t" + pd.covars.get(k - 1));
            }
            validCovars = false;
            invalids.add(pd.covars.get(k - 1));
          } else {
            baseData[k - 1] = d.doubleValue();
            covarData1[k] = d.doubleValue();
          }
        }
        for (int k = 0; k < covarOrder.size(); k++) {
          String covarKey = covarOrder.get(k);
          Double d = covarData.get(covarKey).get(indiv.getKey());
          if (d == null || d.isNaN()) {
            if (writer != null) {
              writer.println(indiv.getKey() + "\t" + covarKey);
            }
            validCovars = false;
            invalids.add(covarKey);
          } else {
            baseData[pd.covars.size() + k] = d.doubleValue();
            covarData1[pd.covars.size() + k + 1] = d.doubleValue();
          }
        }
        if (validCovars) {
          depData.add(pdi.getDepvar());
          baselineIndeps.add(baseData);
          indepData.add(covarData1);
        }
      }
    }
    if (invalids.uniqueSet().size() > 0) {
      for (String k : invalids.uniqueSet()) {
        log.reportWarning("Missing " + invalids.getCount(k) + " for covariate " + k + ".");
      }
    }

    double[][] baseCovars = new double[baselineIndeps.size()][];
    double[][] covars = new double[indepData.size()][];
    for (int k = 0; k < covars.length; k++) {
      covars[k] = indepData.get(k);
      baseCovars[k] = baselineIndeps.get(k);
    }
    int cases = 0;
    int controls = 0;
    Multiset<Double> phenos = ImmutableMultiset.copyOf(depData);
    if (phenos.entrySet().size() == 2) {
      if (phenos.elementSet().equals(ImmutableSet.of(Double.valueOf(0.0), Double.valueOf(1.0)))) {
        cases = phenos.count(Double.valueOf(1.0));
        controls = phenos.count(Double.valueOf(0.0));
      } else if (phenos.elementSet()
                       .equals(ImmutableSet.of(Double.valueOf(1.0), Double.valueOf(2.0)))) {
        cases = phenos.count(Double.valueOf(2.0));
        controls = phenos.count(Double.valueOf(1.0));
      } else {
        log.reportError("Unrecognized case/control designations, cannot count cases and controls");
      }
    }
    RegressionModel baseModel = RegressionModel.determineAppropriate(Doubles.toArray(depData),
                                                                     baseCovars, false, true);
    RegressionModel model = RegressionModel.determineAppropriate(Doubles.toArray(depData), covars,
                                                                 false, true);

    RegressionResult.Builder rr = RegressionResult.builder();
    rr.setNum(depData.size());
    rr.setnCases(cases);
    rr.setnControls(controls);
    if (model.analysisFailed()) {
      rr.setBaseRSq(baseModel.analysisFailed() ? Double.NaN : baseModel.getRsquare());
      rr.setBeta(Double.NaN);
      rr.setSe(Double.NaN);
      rr.setRsq(Double.NaN);
      rr.setPval(Double.NaN);
      rr.setLogistic(model.isLogistic());
    } else {
      int ind = -1;
      for (int l = 0; l < model.getVarNames().length; l++) {
        if ("Indep 1".equals(model.getVarNames()[l])) {
          ind = l;
          break;
        }
      }
      if (ind == -1) {
        rr.setBaseRSq(baseModel.analysisFailed() ? Double.NaN : baseModel.getRsquare());
        rr.setBeta(Double.NaN);
        rr.setSe(Double.NaN);
        rr.setRsq(model.getRsquare());
        rr.setPval(Double.NaN);
        rr.setLogistic(model.isLogistic());
        rr.setStats(Double.NaN);
      } else {
        rr.setBaseRSq(baseModel.analysisFailed() ? Double.NaN : baseModel.getRsquare());
        rr.setBeta(model.getBetas()[ind]);
        rr.setSe(model.getSEofBs()[ind]);
        rr.setRsq(model.getRsquare());
        rr.setPval(model.getSigs()[ind]);
        rr.setLogistic(model.isLogistic());
        rr.setStats(model.getStats()[ind]);
      }
    }
    RegressionResult rrResult = rr.build();
    return rrResult;
  }

  private static class MetaAnalysisResult {
    final double beta;
    final double se;
    final double pval;
    final String direction;
    final int num;
    final int cases;
    final int controls;
    final double i2;
    final double df;
    final double chi;

    public MetaAnalysisResult(double beta, double se, double pval, String direction, int num,
                              int cases, int controls, double i2, double df, double chi) {
      this.beta = beta;
      this.se = se;
      this.pval = pval;
      this.direction = direction;
      this.num = num;
      this.cases = cases;
      this.controls = controls;
      this.i2 = i2;
      this.df = df;
      this.chi = chi;
    }

    @Override
    public String toString() {
      return new StringJoiner("\t").add(Double.toString(beta)).add(Double.toString(se))
                                   .add(Double.toString(pval)).add(direction)
                                   .add(Integer.toString(num)).add(Integer.toString(cases))
                                   .add(Integer.toString(controls)).add(Double.toString(i2))
                                   .add(Double.toString(df)).add(Double.toString(chi)).toString();
    }

  }

  private static class MetaAnalysisSummary {
    MetaAnalysisResult prsResult;
    Map<String, MetaAnalysisResult> perMarkerResultsUni;
    Map<String, MetaAnalysisResult> perMarkerResultsMulti;

    public MetaAnalysisSummary() {
      perMarkerResultsUni = new HashMap<>();
      perMarkerResultsMulti = new HashMap<>();
    }
  }

  private void processMetaAnalyses() {

    if (!Files.exists(metaDir + META_ANALYSIS_DEFINITIONS_FILE)) {
      return;
    }

    // create directory structure
    new File(metaDir + "/" + META_DIRECTORY).mkdir();
    for (MetaFile mf : metaFiles) {
      String dataFolder = metaDir + "/" + META_DIRECTORY + "/" + mf.metaRoot + "/";
      for (Constraint constr : analysisConstraints) {
        String constraintFolder = dataFolder + constr.analysisString + "/";
        for (MetaAnalysis ma : metaAnalyses) {
          String metaFolder = constraintFolder + ma.name + "/";
          File f = new File(metaFolder);
          if (!(f.exists())) {
            f.mkdirs();
          }
        }
      }
    }

    for (MetaAnalysis ma : metaAnalyses) {

      for (MetaFile mf : metaFiles) {

        for (Constraint constr : analysisConstraints) {

          PrintWriter writer = Files.getAppropriateWriter(metaDir + "/" + META_DIRECTORY + "/"
                                                          + ma.name + "_" + mf.metaRoot + "_"
                                                          + constr.analysisString + "_InvVar.out");
          writer.println("SUBTYPE\tVARIATE\tBETA\tSE\tPVAL\tDIRECTION\tNUM\tCASES\tCONTROLS\tHetISq\tHetDF\tHetChiSq");

          final StringBuilder prsDir = new StringBuilder();
          final List<double[]> prsBS = new ArrayList<>();

          final Set<String> allMarkerNamesSet = new HashSet<>();

          ma.entries.forEach(mae -> {
            Study study = mae.study;
            PhenoData pd = mae.study.phenoData.get(mae.phenoFile);

            allMarkerNamesSet.addAll(study.markerRegressions.get(constr, mf).row(pd).keySet());
          });

          final List<String> allMarkerNames = new ArrayList<>(allMarkerNamesSet);
          allMarkerNames.sort(Comparator.naturalOrder());

          int prsCount = 0;
          int prsCase = 0;
          int prsCtrl = 0;
          final Map<String, List<double[]>> perMarkerBS = new HashMap<>();
          final Map<String, List<double[]>> perMarkerMultivarBS = new HashMap<>();
          final Map<String, int[]> perMarkerSums = new HashMap<>();
          final Map<String, int[]> perMarkerMultivarSums = new HashMap<>();
          final Map<String, StringBuilder> perMarkerUniDirections = new HashMap<>();
          final Map<String, StringBuilder> perMarkerMultiDirections = new HashMap<>();

          for (String mkr : allMarkerNames) {
            perMarkerBS.put(mkr, new ArrayList<>());
            perMarkerMultivarBS.put(mkr, new ArrayList<>());
            perMarkerSums.put(mkr, new int[] {0, 0, 0});
            perMarkerMultivarSums.put(mkr, new int[] {0, 0, 0});
            perMarkerUniDirections.put(mkr, new StringBuilder());
            perMarkerMultiDirections.put(mkr, new StringBuilder());
          }

          for (MetaAnalysisEntry mae : ma.entries) {
            Study study = mae.study;
            PhenoData pd = mae.study.phenoData.get(mae.phenoFile);

            RegressionResult rr;
            rr = study.regressions.get(constr, mf).get(pd);
            prsBS.add(new double[] {rr.getBeta(), rr.getSe()});
            prsDir.append(getDirection(rr.getBeta()));
            prsCount += rr.getNum();
            prsCase += rr.getnCases();
            prsCtrl += rr.getnControls();

            Map<String, RegressionResult> mkrRegsUni = study.markerRegressions.get(constr, mf)
                                                                              .row(pd);
            Map<String, RegressionResult> mkrRegsMulti = study.multivarRegressions.get(constr, mf)
                                                                                  .row(pd);

            for (String mkr : allMarkerNames) {
              if (!mkrRegsUni.containsKey(mkr)) {
                perMarkerUniDirections.get(mkr).append("?");
              } else {
                perMarkerUniDirections.get(mkr).append(getDirection(mkrRegsUni.get(mkr).getBeta()));
                perMarkerBS.get(mkr).add(new double[] {mkrRegsUni.get(mkr).getBeta(),
                                                       mkrRegsUni.get(mkr).getSe()});
                perMarkerSums.get(mkr)[0] = perMarkerSums.get(mkr)[0]
                                            + mkrRegsUni.get(mkr).getNum();
                perMarkerSums.get(mkr)[1] = perMarkerSums.get(mkr)[1]
                                            + mkrRegsUni.get(mkr).getnCases();
                perMarkerSums.get(mkr)[2] = perMarkerSums.get(mkr)[2]
                                            + mkrRegsUni.get(mkr).getnControls();
              }
              if (!mkrRegsMulti.containsKey(mkr)) {
                perMarkerMultiDirections.get(mkr).append("?");
              } else {
                perMarkerMultiDirections.get(mkr)
                                        .append(getDirection(mkrRegsMulti.get(mkr).getBeta()));
                perMarkerMultivarBS.get(mkr).add(new double[] {mkrRegsMulti.get(mkr).getBeta(),
                                                               mkrRegsMulti.get(mkr).getSe()});
                perMarkerMultivarSums.get(mkr)[0] = perMarkerMultivarSums.get(mkr)[0]
                                                    + mkrRegsMulti.get(mkr).getNum();
                perMarkerMultivarSums.get(mkr)[1] = perMarkerMultivarSums.get(mkr)[1]
                                                    + mkrRegsMulti.get(mkr).getnCases();
                perMarkerMultivarSums.get(mkr)[2] = perMarkerMultivarSums.get(mkr)[2]
                                                    + mkrRegsMulti.get(mkr).getnControls();
              }
            }
          }

          MetaAnalysisSummary summary = new MetaAnalysisSummary();

          // Meta-analyze each snp across studies
          double[] data = org.pankratzlab.utils.gwas.MetaAnalysis.inverseVarianceWeighting(prsBS.toArray(new double[0][]),
                                                                                           log);
          double pval = computePVal(data[0], data[1], prsCount);

          int df = ma.entries.size() - 1;
          Map<String, Double> perMarkerUniI2 = new HashMap<>();
          Map<String, Double> perMarkerMultiI2 = new HashMap<>();
          Map<String, Double> perMarkerHetStat = new HashMap<>();
          Map<String, Double> perMarkerMultiHetStat = new HashMap<>();
          perMarkerBS.entrySet().stream().forEach(betaSEs -> {
            String marker = betaSEs.getKey();
            double[] hetStati2 = computeHetStat(df, betaSEs.getValue());
            perMarkerHetStat.put(marker, hetStati2[0]);
            perMarkerUniI2.put(marker, hetStati2[1]);
          });
          perMarkerMultivarBS.entrySet().stream().forEach(betaSEs -> {
            String marker = betaSEs.getKey();
            double[] hetStati2 = computeHetStat(df, betaSEs.getValue());
            perMarkerMultiHetStat.put(marker, hetStati2[0]);
            perMarkerMultiI2.put(marker, hetStati2[1]);
          });

          summary.prsResult = new MetaAnalysisResult(data[0], data[1], pval, prsDir.toString(),
                                                     prsCount, prsCase, prsCtrl, Double.NaN,
                                                     Double.NaN, Double.NaN);

          writer.println("PRS\tPRS\t" + summary.prsResult.toString());

          for (Entry<String, List<double[]>> cell : perMarkerBS.entrySet()) {
            double[] mkrData = org.pankratzlab.utils.gwas.MetaAnalysis.inverseVarianceWeighting(cell.getValue()
                                                                                                    .toArray(new double[0][]),
                                                                                                log);
            double mkrPval = computePVal(mkrData[0], mkrData[1],
                                         perMarkerSums.get(cell.getKey())[0]);
            double chi = ProbDist.ChiDist(perMarkerHetStat.get(cell.getKey()), df);
            MetaAnalysisResult mar = new MetaAnalysisResult(mkrData[0], mkrData[1], mkrPval,
                                                            perMarkerUniDirections.get(cell.getKey())
                                                                                  .toString(),
                                                            perMarkerSums.get(cell.getKey())[0],
                                                            perMarkerSums.get(cell.getKey())[1],
                                                            perMarkerSums.get(cell.getKey())[2],
                                                            perMarkerUniI2.get(cell.getKey()), df,
                                                            chi);
            summary.perMarkerResultsUni.put(cell.getKey(), mar);
            writer.println(cell.getKey() + "\tUNIVARIATE\t" + mar.toString());
          }
          for (Entry<String, List<double[]>> cell : perMarkerMultivarBS.entrySet()) {
            double[] mkrData = org.pankratzlab.utils.gwas.MetaAnalysis.inverseVarianceWeighting(cell.getValue()
                                                                                                    .toArray(new double[0][]),
                                                                                                log);
            double mkrPval = computePVal(mkrData[0], mkrData[1],
                                         perMarkerMultivarSums.get(cell.getKey())[0]);
            double chi = ProbDist.ChiDist(perMarkerMultiHetStat.get(cell.getKey()), df);
            MetaAnalysisResult mar = new MetaAnalysisResult(mkrData[0], mkrData[1], mkrPval,
                                                            perMarkerMultiDirections.get(cell.getKey())
                                                                                    .toString(),
                                                            perMarkerMultivarSums.get(cell.getKey())[0],
                                                            perMarkerMultivarSums.get(cell.getKey())[1],
                                                            perMarkerMultivarSums.get(cell.getKey())[2],
                                                            perMarkerMultiI2.get(cell.getKey()), df,
                                                            chi);
            summary.perMarkerResultsMulti.put(cell.getKey(), mar);
            writer.println(cell.getKey() + "\tMULTIVARIATE\t" + mar.toString());
          }

          writer.close();

          summaryResults.get(constr, mf).put(ma, summary);
        }
      }
    }
  }

  private double[] computeHetStat(int df, List<double[]> betaSes) {
    double[] z = new double[betaSes.size()];
    double[] w = new double[betaSes.size()];
    double stat = 0;
    double weight = 0;

    for (int i = 0; i < betaSes.size(); i++) {
      z[i] = betaSes.get(i)[0];
      w[i] = 1d / (betaSes.get(i)[1] * betaSes.get(i)[1]);
      stat += (w[i] * z[i]);
      weight += w[i];
    }

    double e = stat / weight;

    double hetStat = 0;
    for (int i = 0; i < z.length; i++) {
      hetStat += (z[i] - e) * (z[i] - e) * w[i];
    }

    double i2 = (hetStat - df + 1) / (hetStat * 100d);
    return new double[] {hetStat, i2};
  }

  private double computePVal(double beta, double se, double n) {
    return ProbDist.NormDist(Math.abs(beta / se));
  }

  private String getExcelPVal(double beta, double se) {
    return "=(1-(NORM.S.DIST(ABS(" + beta + "/" + se + "),TRUE)))*2";
  }

  private void writeResults() {
    String resFile = metaDir + "results.xln";
    log.report("Writing regression results... [ --> " + resFile + "]");
    try (PrintWriter writer = Files.getAppropriateWriter(resFile)) {
      writer.println(REGRESSION_HEADER);

      for (Study study : studies) {
        for (MetaFile mf : metaFiles) {
          String dFile = mf.metaFile;
          String dataFile = mf.metaRoot;
          for (Constraint constr : analysisConstraints) {
            String key = constr.analysisString;
            TrioScoreTest.Results trioTestResults = study.trioScoreTests.get(constr, mf);
            final String pairedStatTrios;
            final String pairedT;
            final String wilcoxon;
            if (trioTestResults == null) {
              pairedStatTrios = "0";
              pairedT = "N/A";
              wilcoxon = "N/A";
            } else {
              pairedStatTrios = String.valueOf(trioTestResults.getTrioCount());
              pairedT = String.valueOf(trioTestResults.pairedTPVal);
              wilcoxon = String.valueOf(trioTestResults.wilcoxonPVal);
            }

            String resultPrefix = new StringJoiner("\t").add(study.studyName).add(dataFile)
                                                        .add(ext.formSciNot(constr.indexThreshold,
                                                                            5, false))
                                                        .toString();

            String middle = new StringJoiner("\t").add(pairedT).add(wilcoxon).add(pairedStatTrios)
                                                  .toString();
            String endPRS = new StringJoiner("\t").add(String.valueOf(dataCounts.get(dFile)
                                                                                .get(key)))
                                                  .add(String.valueOf(study.hitWindowCnts.get(constr,
                                                                                              mf)))
                                                  .add(String.valueOf(study.hitSnpCounts.get(constr,
                                                                                             mf)))
                                                  .add(String.valueOf(study.scores.get(constr,
                                                                                       mf)[0]))
                                                  .add(String.valueOf(study.scores.get(constr,
                                                                                       mf)[1]))
                                                  .toString();
            String endUNIVAR = new StringJoiner("\t").add("1").add("1").add("1")
                                                     .add(String.valueOf(study.scores.get(constr,
                                                                                          mf)[0]))
                                                     .add(String.valueOf(study.scores.get(constr,
                                                                                          mf)[1]))
                                                     .toString();

            if (study.phenoFiles.isEmpty()) {
              writeSingleResult("--", "PRS", "UNIVARIATE", RegressionResult.dummy(), resultPrefix,
                                middle, endPRS, writer);
            } else {
              for (String pheno : study.phenoFiles) {
                RegressionResult rr = study.regressions.get(constr, mf)
                                                       .get(study.phenoData.get(pheno));
                if (rr == null) {
                  rr = RegressionResult.dummy();
                }
                String phenoName = pheno.endsWith(".pheno") ? pheno.substring(0, pheno.length() - 6)
                                                            : pheno;

                writeSingleResult(phenoName, "PRS", "UNIVARIATE", rr, resultPrefix, middle, endPRS,
                                  writer);

                Set<String> mkrs = new LinkedHashSet<>();
                mkrs.addAll(study.markerRegressions.get(constr, mf).columnKeySet());
                mkrs.addAll(study.multivarRegressions.get(constr, mf).columnKeySet());

                Map<String, RegressionResult> res = study.markerRegressions.get(constr, mf)
                                                                           .row(study.phenoData.get(pheno));
                mkrs.stream().filter(res::containsKey).forEach(s -> {
                  writeSingleResult(phenoName, s, "UNIVARIATE", res.get(s), resultPrefix, middle,
                                    endUNIVAR, writer);
                });
                Map<String, RegressionResult> resMult = study.multivarRegressions.get(constr, mf)
                                                                                 .row(study.phenoData.get(pheno));
                mkrs.stream().filter(resMult::containsKey).forEach(s -> {
                  final RegressionResult rr2 = resMult.get(s);
                  String myEndPRS = new StringJoiner("\t").add(String.valueOf(dataCounts.get(dFile)
                                                                                        .get(key)
                                                                              - rr2.droppedSnps))
                                                          .add(String.valueOf(study.hitWindowCnts.get(constr,
                                                                                                      mf)
                                                                              - rr2.droppedSnps))
                                                          .add(String.valueOf(study.hitSnpCounts.get(constr,
                                                                                                     mf)))
                                                          .add(String.valueOf(study.scores.get(constr,
                                                                                               mf)[0]))
                                                          .add(String.valueOf(study.scores.get(constr,
                                                                                               mf)[1]))
                                                          .toString();
                  writeSingleResult(phenoName, s, "MULTIVARIATE", rr2, resultPrefix, middle,
                                    myEndPRS, writer);
                });
              }
            }
          }
        }
      }

      for (MetaAnalysis ma : metaAnalyses) {
        for (MetaFile mf : metaFiles) {
          String dataFile = mf.metaRoot;
          for (Constraint constr : analysisConstraints) {
            MetaAnalysisSummary summary = summaryResults.get(constr, mf).get(ma);

            final MetaAnalysisResult result = summary.prsResult;

            String linStr = formResultLine(ma, dataFile, constr, result, "PRS", true);

            writer.println(linStr);

            for (Entry<String, MetaAnalysisResult> res : summary.perMarkerResultsUni.entrySet()) {
              writer.println(formResultLine(ma, dataFile, constr, res.getValue(), res.getKey(),
                                            true));
            }

            for (Entry<String, MetaAnalysisResult> res : summary.perMarkerResultsMulti.entrySet()) {
              writer.println(formResultLine(ma, dataFile, constr, res.getValue(), res.getKey(),
                                            false));
            }

          }
        }
      }
    }
  }

  private String formResultLine(MetaAnalysis ma, String dataFile, Constraint constr,
                                final MetaAnalysisResult result, String name, boolean univar) {
    StringJoiner line = new StringJoiner("\t");
    line.add("META");
    line.add(dataFile);
    line.add(ext.formSciNot(constr.indexThreshold, 5, false));
    line.add(ma.name);
    line.add(name);
    line.add(univar ? "UNIVARIATE" : "MULTIVARIATE");
    line.add(".");
    line.add(".");
    line.add(".");
    line.add(Double.toString(result.pval));
    line.add(getExcelPVal(result.beta, result.se));
    line.add(Double.toString(result.beta));
    line.add(Double.toString(result.se));
    line.add(result.direction);
    line.add(Integer.toString(result.num));
    line.add(Integer.toString(result.cases));
    line.add(Integer.toString(result.controls));
    for (int i = 0; i < 8; i++) {
      line.add(".");
    }
    String linStr = line.toString();
    return linStr;
  }

  private void writeAndRunMR() {
    Map<MR_Key, String> mrrScripts = new HashMap<>();
    Map<MR_Key, String> metaMrrScripts = new HashMap<>();
    for (Study study : studies) {
      mrrScripts.putAll(writeMarkerResults(study));
    }
    if (!metaAnalyses.isEmpty()) {
      metaMrrScripts.putAll(writeMetaMarkerResults());
    }
    for (String script : mrrScripts.values()) {
      log.report("Executing " + script);
      CmdLine.basic(log).run(Command.basic(script));
    }
    if (!metaAnalyses.isEmpty()) {
      log.report("Executing meta-analysis MR scripts...");
      for (String script : metaMrrScripts.values()) {
        log.report("Executing " + script);
        CmdLine.basic(log).run(Command.basic(script));
      }
    }

    parseMRResults(mrrScripts, metaMrrScripts);

  }

  private void parseMRResults(Map<MR_Key, String> mrrScripts, Map<MR_Key, String> metaMrrScripts) {
    AtomicBoolean header = new AtomicBoolean(true);
    AtomicBoolean append = new AtomicBoolean(false);
    Streams.concat(mrrScripts.keySet().stream(), metaMrrScripts.keySet().stream()).forEach(key -> {
      parseMRResultFile(header, append, key);
    });
  }

  private void parseMRResultFile(AtomicBoolean header, AtomicBoolean append, MR_Key key) {
    String dir;
    if (key.isMeta()) {
      dir = metaDir + META_DIRECTORY + "/" + key.mf.metaRoot + "/" + key.c.analysisString + "/"
            + key.name + "/" + key.subdir;
    } else {
      dir = getDirPath(key.study, key.mf, key.c) + key.subdir + "/" + key.name;
    }
    String[] resFiles = Files.list(dir, key.subtype.resultFile);
    if (resFiles.length == 0) {
      log.reportWarning("No MendelianRandomization results found in " + dir);
    } else if (resFiles.length > 1) {
      log.reportWarning("Multiple MendelianRandomization results found in " + dir);
    } else {
      String resFile = resFiles[0];

      FileColumn<String> studyCol = new FixedValueColumn("STUDY",
                                                         key.isMeta() ? "META"
                                                                      : key.study.studyName);
      FileColumn<String> fileCol = new FixedValueColumn("DATAFILE", key.mf.metaRoot);
      FileColumn<String> phenoCol = new FixedValueColumn("PHENOTYPE", key.name);
      FileColumn<String> constrCol = new FixedValueColumn("CONSTRAINT", key.c.analysisString);
      FileColumn<String> subtypeCol = new FixedValueColumn("SUBTYPE", key.subtype.name());
      FileColumn<String> variateCol = new FixedValueColumn("VARIATE", key.multi ? "MULTIVARIATE"
                                                                                : "UNIVARIATE");

      FileColumn<String> methodCol = new AliasedFileColumn("MR METHOD", key.subtype.methodCol);
      FileColumn<Double> betaCol = new DoubleWrapperColumn("BETA",
                                                           new AliasedFileColumn(key.subtype.betaCol));
      FileColumn<Double> seCol = new DoubleWrapperColumn("SE",
                                                         new AliasedFileColumn(key.subtype.seCol));
      FileColumn<Double> pCol = new DoubleWrapperColumn("P-VALUE",
                                                        new AliasedFileColumn(key.subtype.pCol));
      FileColumn<?> nCol = key.subtype.hasNCol() ? new IntegerWrapperColumn("NUM",
                                                                            new AliasedFileColumn(key.subtype.nCol))
                                                 : new FixedValueColumn("NUM", ".");

      String sigInMeta = String.valueOf(dataCounts.get(key.mf.metaFile).get(key.c.analysisString));
      String indexVarMeta = key.isMeta() ? "."
                                         : (key.multi ? String.valueOf(key.study.hitWindowCnts.get(key.c,
                                                                                                   key.mf))
                                                      : "1");
      String indexVarData = key.isMeta() ? "."
                                         : (key.multi ? String.valueOf(key.study.hitSnpCounts.get(key.c,
                                                                                                  key.mf))
                                                      : "1");

      FileColumn<String> nCol2 = new FixedValueColumn("#sigInMeta", sigInMeta);
      FileColumn<String> nCol3 = new FixedValueColumn("#indexVariantsInMeta", indexVarMeta);
      FileColumn<String> nCol4 = new FixedValueColumn("#indexVariantsInDataset", indexVarData);

      List<FileColumn<?>> cols = Lists.of(studyCol, fileCol, phenoCol, constrCol, subtypeCol,
                                          variateCol, methodCol, betaCol, seCol, pCol, nCol, nCol2,
                                          nCol3, nCol4);

      try {
        FileParserFactory.setup(dir + "/" + resFile, cols.toArray(FileColumn[]::new)).build()
                         .parseToFile(GeneScorePipeline.this.metaDir + "mr_results.xln", "\t",
                                      header.getAndSet(false), append.getAndSet(true));
      } catch (IOException e1) {
        throw new RuntimeException(e1);
      }

    }
  }

  private void writeSingleResult(String pheno, String subtype, String variate, RegressionResult rr,
                                 String resultPrefix, String middle, String end,
                                 PrintWriter writer) {
    String pvalExcl = rr.getNum() == 0 ? "."
                                       : (rr.isLogistic() ? getExcelPVal(rr.getBeta(), rr.getSe())
                                                          : "=TDIST(" + Math.abs(rr.getStats())
                                                            + "," + rr.getNum() + ",2)");

    String line = Joiner.on("\t")
                        .join(ImmutableList.builder().add(resultPrefix).add(pheno).add(subtype)
                                           .add(variate).add(rr.getBaseRSq()).add(rr.getRsq())
                                           .add((Double.isNaN(rr.getRsq()) ? Double.NaN
                                                                           : (Double.isNaN(rr.getBaseRSq()) ? rr.getRsq()
                                                                                                            : (new BigDecimal(rr.getRsq()
                                                                                                                              + "")).subtract(new BigDecimal(rr.getBaseRSq()
                                                                                                                                                             + "")))))
                                           .add(rr.getPval()).add(pvalExcl).add(rr.getBeta())
                                           .add(rr.getSe()).add(getDirection(rr.getBeta()))
                                           .add(rr.getNum()).add(rr.getnCases())
                                           .add(rr.getnControls()).add(middle).add(end).build());

    writer.println(line);

  }

  private String getDirection(double beta) {
    if (Double.isNaN(beta)) return "?";
    if (beta > 0) return "+";
    if (beta < 0) return "-";
    if (beta == 0) return "0";
    return "?";
  }

  enum MR_TYPE {
    MR(".csv", "Method", "Estimate", "Std Error", "P-value", null),
    TWO_SAMP("MR_results.txt", "method", "b", "se", "pval", "nsnp");

    final String resultFile;
    final String methodCol;
    final String betaCol;
    final String seCol;
    final String pCol;
    final String nCol;

    private MR_TYPE(String resultFile, String methodCol, String betaCol, String seCol, String pCol,
                    String nCol) {
      this.resultFile = resultFile;
      this.methodCol = methodCol;
      this.betaCol = betaCol;
      this.seCol = seCol;
      this.pCol = pCol;
      this.nCol = nCol;
    }

    boolean hasNCol() {
      return nCol != null;
    }

  }

  static class MR_Key {
    Study study;
    MetaFile mf;
    Constraint c;
    String name;
    boolean multi;
    MR_TYPE subtype;
    String subdir;

    public MR_Key(Study study, MetaFile mf, Constraint c, String name, boolean multi,
                  MR_TYPE subtype, String subdir) {
      this.study = study;
      this.mf = mf;
      this.c = c;
      this.name = name;
      this.multi = multi;
      this.subtype = subtype;
      this.subdir = subdir;
    }

    boolean isMeta() {
      return study == null;
    }

  }

  private Map<MR_Key, String> writeMetaMarkerResults() {
    Map<MR_Key, String> mrrScripts = new HashMap<>();

    for (MetaFile mf : metaFiles) {
      String dataFile = mf.metaRoot;
      for (Constraint constr : analysisConstraints) {
        String analysisKey = constr.analysisString;
        for (MetaAnalysis ma : metaAnalyses) {
          String prefDir = metaDir + "/" + META_DIRECTORY + "/" + dataFile + "/" + analysisKey + "/"
                           + ma.name + "/";
          writeMetaRScripts(mrrScripts, mf, dataFile, constr, ma, prefDir, true);
          writeMetaRScripts(mrrScripts, mf, dataFile, constr, ma, prefDir, false);
        }

      }
    }

    return mrrScripts;
  }

  private void writeMetaRScripts(Map<MR_Key, String> mrrScripts, MetaFile mf, String dataFile,
                                 Constraint constr, MetaAnalysis ma, String prefDir,
                                 boolean univar) {
    MetaAnalysisSummary summary = summaryResults.get(constr, mf).get(ma);

    String mrPrefDirMRStr = prefDir + "/mr_" + (univar ? "univar" : "multivar") + "/";
    File mrPrefDirMR = new File(mrPrefDirMRStr);
    mrPrefDirMR.mkdirs();
    String mrPrefDirTSStr = prefDir + "/two-sample_" + (univar ? "univar" : "multivar") + "/";
    File mrPrefDirTS = new File(mrPrefDirTSStr);
    mrPrefDirTS.mkdirs();

    PrintWriter twoSampleExposure = Files.getAppropriateWriter(mrPrefDirTSStr + "exposure_"
                                                               + ma.name + ".out");
    PrintWriter twoSampleOutcome = Files.getAppropriateWriter(mrPrefDirTSStr + "outcome_" + ma.name
                                                              + ".out");
    twoSampleExposure.println(TWO_SAMPLE_HDR);
    twoSampleOutcome.println(TWO_SAMPLE_HDR);

    // exposure is info from meta file
    // outcome is from gsp

    String mrrInputFile = MARKER_REGRESSION_PREFIX + ma.name + MARKER_REGRESSION_EXTEN;
    PrintWriter markerWriter = Files.getAppropriateWriter(mrPrefDirMRStr + mrrInputFile);
    markerWriter.println(MARKER_RESULT_HEADER);

    String resultPrefix = new StringJoiner("\t").add("META").add(dataFile)
                                                .add(ext.formSciNot(constr.indexThreshold, 5,
                                                                    false))
                                                .toString();

    int foundMkrs = 0;
    for (String marker : mf.metaMarkers.keySet()) {
      final Map<String, MetaAnalysisResult> perMarkerResults = univar ? summary.perMarkerResultsUni
                                                                      : summary.perMarkerResultsMulti;
      if (!perMarkerResults.containsKey(marker)) {
        continue;
      }

      foundMkrs++;

      MetaMarker metaMarker = mf.metaMarkers.get(marker);
      double metaBeta = metaMarker.beta;
      double metaSE = metaMarker.se;

      double markerBeta = perMarkerResults.get(marker).beta;
      double markerSE = perMarkerResults.get(marker).se;

      if (!Double.isNaN(markerBeta) && !Double.isNaN(markerSE)) {
        markerWriter.println(resultPrefix + "\t" + ma.name + "\t" + marker + "\t" + metaBeta + "\t"
                             + metaSE + "\t" + markerBeta + "\t" + markerSE);

        twoSampleExposure.println(marker + "\t" + metaMarker.alleles.a1 + "\t"
                                  + metaMarker.alleles.a2 + "\t" + metaMarker.freq + "\t" + metaBeta
                                  + "\t" + metaSE + "\t" + metaMarker.pval + "\tExposure");

        String FREQ = ".";
        twoSampleOutcome.println(marker + "\t" + metaMarker.alleles.a1 + "\t"
                                 + metaMarker.alleles.a2 + "\t" + FREQ + "\t" + markerBeta + "\t"
                                 + markerSE + "\t" + perMarkerResults.get(marker).pval
                                 + "\tOutcome");
      }
    }
    twoSampleExposure.close();
    twoSampleOutcome.close();
    markerWriter.close();
    if (foundMkrs > 1) {
      MR_Key mrKey = new MR_Key(null, mf, constr, ma.name, !univar, MR_TYPE.MR,
                                "mr_" + (univar ? "univar" : "multivar"));
      String mrrScript = writeMRRScript(mrPrefDirMR, ma.name);
      mrrScripts.put(mrKey, mrrScript);

      MR_Key twoSampKey = new MR_Key(null, mf, constr, ma.name, !univar, MR_TYPE.TWO_SAMP,
                                     "two-sample_" + (univar ? "univar" : "multivar"));
      String tsRScript = writeTwoSampleRScript(mrPrefDirTS, ma.name);
      mrrScripts.put(twoSampKey, tsRScript);
    }
  }

  private Map<MR_Key, String> writeMarkerResults(Study study) {
    Map<MR_Key, String> mrrScripts = new HashMap<>();

    for (MetaFile mf : metaFiles) {
      String dataFile = mf.metaRoot;
      for (Constraint constr : analysisConstraints) {
        String analysisKey = constr.analysisString;
        File prefDir = new File(getDirPath(study, mf, constr));
        if (!prefDir.exists()) {
          log.report("Error - no subfolder for '" + analysisKey + "' analysis");
          continue;
        }
        log.report("Writing marker-specific results for "
                   + new StringJoiner("//").add(study.studyName).add(dataFile)
                                           .add(ext.formSciNot(constr.indexThreshold, 5, false))
                                           .toString());
        String resultPrefix = new StringJoiner("\t").add(study.studyName).add(dataFile)
                                                    .add(ext.formSciNot(constr.indexThreshold, 5,
                                                                        false))
                                                    .toString();

        Set<String> mkrs = study.matchedMarkerIndices.get(constr, mf).keySet();
        Set<String> intersectMkrs = Sets.intersection(mkrs, mf.metaMarkers.keySet());
        List<String> markersInOrder = new ArrayList<>(intersectMkrs);
        final Map<String, Integer> map = study.matchedMarkerIndices.get(constr, mf);
        markersInOrder.sort((m1, m2) -> {
          return map.getOrDefault(m1, -1).compareTo(map.getOrDefault(m2, -1));
        });

        for (int i = 0; i < study.phenoFiles.size(); i++) {
          String pheno = study.phenoFiles.get(i);
          writeRScripts(study, mrrScripts, mf, constr, resultPrefix, markersInOrder, pheno, true);
          writeRScripts(study, mrrScripts, mf, constr, resultPrefix, markersInOrder, pheno, false);
        }

        PrintWriter forestWriter = Files.getAppropriateWriter(prefDir + "/forest_input.txt");
        PrintWriter betaWriter = Files.getAppropriateWriter(prefDir + "/" + REGRESSION_BETA_FILE);
        betaWriter.print("MarkerName\tbeta\tse");
        for (String marker : markersInOrder) {
          betaWriter.print("\tbeta.");
          betaWriter.print(marker);
          betaWriter.print("\tse.");
          betaWriter.print(marker);
        }
        betaWriter.println();
        for (int i = 0; i < study.phenoFiles.size(); i++) {
          String pheno = study.phenoFiles.get(i);
          PhenoData pd = study.phenoData.get(pheno);
          RegressionResult rrPheno = study.regressions.get(constr, mf).get(pd);
          betaWriter.print(pd.phenoName);
          betaWriter.print("\t");
          betaWriter.print(rrPheno.getBeta());
          betaWriter.print("\t");
          betaWriter.print(rrPheno.getSe());

          String[][] ids = study.data.get(constr, mf).getIds();

          for (String marker : markersInOrder) {
            int mkrInd = study.matchedMarkerIndices.get(constr, mf).getOrDefault(marker, -1);
            if (mkrInd != -1) {
              Map<String, Double> dosages = new HashMap<>();
              for (int s = 0; s < ids.length; s++) {
                dosages.put(ids[s][0] + "\t" + ids[s][1],
                            study.matchedMarkerDosages.get(constr, mf)[s][mkrInd]);
              }
              RegressionResult rrResult = actualRegression(dosages, null, pd);
              dosages = null;

              betaWriter.print("\t");
              betaWriter.print(rrResult.getBeta());
              betaWriter.print("\t");
              betaWriter.print(rrResult.getSe());
            } else {
              betaWriter.print("\t");
              betaWriter.print(".");
              betaWriter.print("\t");
              betaWriter.print(".");
            }
          }
          betaWriter.println();

          forestWriter.println(pd.phenoName + "\t" + prefDir + "/" + REGRESSION_BETA_FILE + "\t\t"
                               + ext.rootOf(pheno));
        }
        betaWriter.close();
        forestWriter.close();

        new ForestPlot(Optional.of(prefDir + "/forest_input.txt"), log).screenCapAll("plots/",
                                                                                     plotOddsRatio,
                                                                                     true);
      }
    }
    return mrrScripts;
  }

  private void writeRScripts(Study study, Map<MR_Key, String> mrrScripts, MetaFile mf,
                             Constraint constr, String resultPrefix, List<String> markersInOrder,
                             String pheno, boolean univariate) {
    PhenoData pd = study.phenoData.get(pheno);

    String var = univariate ? "univar" : "multivar";
    String subdir_mr = "mr/mr_" + var;
    String subdir_2s = "mr/two-sample_" + var;

    String mrPrefDirMRStr = getDirPath(study, mf, constr) + subdir_mr + "/" + pheno + "/";
    File mrPrefDirMR = new File(mrPrefDirMRStr);
    mrPrefDirMR.mkdirs();
    String mrPrefDirTSStr = getDirPath(study, mf, constr) + subdir_2s + "/" + pheno + "/";
    File mrPrefDirTS = new File(mrPrefDirTSStr);
    mrPrefDirTS.mkdirs();

    PrintWriter twoSampleExposure = Files.getAppropriateWriter(mrPrefDirTSStr + "exposure_" + pheno
                                                               + ".out");
    PrintWriter twoSampleOutcome = Files.getAppropriateWriter(mrPrefDirTSStr + "outcome_" + pheno
                                                              + ".out");
    twoSampleExposure.println(TWO_SAMPLE_HDR);
    twoSampleOutcome.println(TWO_SAMPLE_HDR);

    // exposure is info from meta file
    // outcome is from gsp

    String mrrInputFile = MARKER_REGRESSION_PREFIX + pheno + MARKER_REGRESSION_EXTEN;
    PrintWriter markerWriter = Files.getAppropriateWriter(mrPrefDirMRStr + mrrInputFile);
    markerWriter.println(MARKER_RESULT_HEADER);

    Set<String> droppedNoRegress = new HashSet<>();
    Set<String> droppedNaN = new HashSet<>();
    Set<String> droppedMAC = new HashSet<>();

    int writtenMarkers = 0;
    for (String marker : markersInOrder) {
      if (mf.metaMarkers.containsKey(marker)) { // may have been dropped by meta HitWindows
        RegressionResult rrResult = (univariate ? study.markerRegressions
                                                : study.multivarRegressions).get(constr, mf)
                                                                            .get(pd, marker);

        if (rrResult == null) {
          droppedNoRegress.add(marker);
          continue;
        }

        if (!study.mkrPassesMAC(constr, mf, pd, marker)) {
          droppedMAC.add(marker);
          // filtered by cmac
          continue;
        }

        MetaMarker metaMarker = mf.metaMarkers.get(marker);
        double metaBeta = metaMarker.beta;
        double metaSE = metaMarker.se;

        double markerBeta = rrResult.getBeta();
        double markerSE = rrResult.se;

        if (!Double.isNaN(markerBeta) && !Double.isNaN(markerSE)) {
          writtenMarkers++;
          markerWriter.println(resultPrefix + "\t" + pheno + "\t" + marker + "\t" + metaBeta + "\t"
                               + metaSE + "\t" + markerBeta + "\t" + markerSE);

          twoSampleExposure.println(marker + "\t" + metaMarker.alleles.a1 + "\t"
                                    + metaMarker.alleles.a2 + "\t" + metaMarker.freq + "\t"
                                    + metaBeta + "\t" + metaSE + "\t" + metaMarker.pval
                                    + "\tExposure");

          String FREQ = ".";
          twoSampleOutcome.println(marker + "\t" + metaMarker.alleles.a1 + "\t"
                                   + metaMarker.alleles.a2 + "\t" + FREQ + "\t" + markerBeta + "\t"
                                   + markerSE + "\t" + rrResult.getPval() + "\tOutcome");
        } else {
          droppedNaN.add(marker);
        }
      }
    }
    twoSampleExposure.close();
    twoSampleOutcome.close();
    markerWriter.close();
    if (writtenMarkers > 0) {
      MR_Key mrKey = new MR_Key(study, mf, constr, pheno, !univariate, MR_TYPE.MR, subdir_mr);
      String mrrScript = writeMRRScript(mrPrefDirMR, pheno);
      mrrScripts.put(mrKey, mrrScript);

      MR_Key twoSKey = new MR_Key(study, mf, constr, pheno, !univariate, MR_TYPE.TWO_SAMP,
                                  subdir_2s);
      String tsRScript = writeTwoSampleRScript(mrPrefDirTS, pheno);
      mrrScripts.put(twoSKey, tsRScript);
    }

    if (droppedNaN.size() > 0) {
      log.reportWarning("Dropped " + droppedNaN.size() + " markers ["
                        + ArrayUtils.toStr(droppedNaN, ", ")
                        + "] for NaN regression results (beta/se) when writing MR input files for "
                        + mrPrefDirMRStr);
    }
    if (droppedMAC.size() > 0) {
      log.reportWarning("Dropped " + droppedMAC.size() + " markers ["
                        + ArrayUtils.toStr(droppedMAC, ", ")
                        + "] for not passing MAC filter when writing MR input files for "
                        + mrPrefDirMRStr);
    }
    if (droppedNoRegress.size() > 0) {
      log.reportWarning("Dropped " + droppedNoRegress.size() + " markers ["
                        + ArrayUtils.toStr(droppedNoRegress, ", ")
                        + "] for missing regression results when writing MR input files for "
                        + mrPrefDirMRStr);
    }
  }

  private String writeTwoSampleRScript(File prefDir, String pheno) {
    List<String> commands = new ArrayList<>();
    commands.add("options(repos = c(CRAN = \"https://cloud.r-project.org\"))");
    commands.add("setwd(\"" + ext.verifyDirFormat(prefDir.getAbsolutePath()) + "\")");
    commands.add("if (!require(devtools)) {");
    commands.add("  install.packages(\"devtools\")");
    commands.add("}");
    commands.add("if (!require(markdown)) {");
    commands.add("  install.packages(\"markdown\")");
    commands.add("}");
    commands.add("if (!require(TwoSampleMR)) {");
    commands.add("  devtools::install_github(\"MRCIEU/TwoSampleMR\")");
    commands.add("}");
    commands.add("library(TwoSampleMR)");

    commands.add("exp <- read_exposure_data(filename=\"exposure_" + pheno + ".out\", sep='\\t')");
    commands.add("out <- read_outcome_data(filename=\"outcome_" + pheno + ".out\", sep='\\t')");
    commands.add("dat <-harmonise_data(exp, out, action = 1)");
    commands.add("mr_results <-mr(dat)");
    commands.add("write.table(mr_results, \"MR_results.txt\", sep=\"\t\", col.names=T, row.names=F, quote=F)");
    commands.add("mr_report(dat)");
    commands.add("res_single <- mr_singlesnp(dat, all_method = c(\"mr_weighted_median\", \"mr_egger_regression\", \"mr_ivw\"))");
    commands.add("mr_forest <- mr_forest_plot(res_single)");
    commands.add("png('mr_forest_plot.png')");
    commands.add("mr_forest[[1]]");
    commands.add("dev.off()");

    commands.add("res_single <- mr_singlesnp(dat, all_method = c(\"mr_ivw\", \"mr_egger_regression\", \"mr_weighted_median\"))");
    commands.add("write.table(res_single, \"MRsingleSNP_" + pheno
                 + ".txt\", sep=\"\\t\", col.names=T, row.names=F, quote=F)");

    commands.add("loo <- mr_leaveoneout(dat)");
    commands.add("write.table(loo, \"MR_loo_" + pheno
                 + ".txt\", sep=\"\\t\", col.names=T, row.names=F, quote=F)");

    commands.add("egger.radial<-RadialMR::egger_radial(dat)");
    commands.add("egger.outliers<- egger.radial$outliers");
    commands.add("write.table(egger.outliers, \"MR_egger_outlier_" + pheno
                 + ".txt\", sep=\"\\t\", col.names=T, row.names=F, quote=F)");

    // commands.add("ivw.radial<-RadialMR::ivw_radial(dat)");
    // commands.add("ivw.outliers<- ivw.radial$outliers");
    // commands.add("write.table(ivw.outliers, \"MR_ivw_outlier_" + pheno
    // + ".txt\", sep=\"\\t\", col.names=T, row.names=F, quote=F)");
    //
    // commands.add("png(filename='RadialMR.png')");
    // commands.add("RadialMR::plot_radial(c(ivw.radial, egger.radial),T,F,F)");
    // commands.add("dev.off()");

    Files.writeIterable(commands, ext.verifyDirFormat(prefDir.getAbsolutePath()) + TWO_SAMPLE_PREFIX
                                  + pheno + ".R");

    String setEnvCmd = Files.isWindows() ? "set" : "export";
    commands = new ArrayList<>();
    commands.add("cd " + prefDir.getAbsolutePath());
    if (rLibsDir != null) {
      commands.add(setEnvCmd + " R_LIBS=" + rLibsDir);
      commands.add(setEnvCmd + " R_LIBS_SITE=" + rLibsDir);
    } else {
      log.report("No R library directory specified, MendelianRandomization library will be installed if not present in default R library directory.");
    }
    commands.add("Rscript " + ext.verifyDirFormat(prefDir.getAbsolutePath()) + TWO_SAMPLE_PREFIX
                 + pheno + ".R");

    String runScript = ext.verifyDirFormat(prefDir.getAbsolutePath()) + TWO_SAMPLE_PREFIX + pheno
                       + (Files.isWindows() ? ".bat" : ".sh");
    Files.writeIterable(commands, runScript);
    Files.chmod(runScript);
    return runScript;
  }

  private String writeMRRScript(File prefDir, String pheno) {

    List<String> commands = new ArrayList<>();
    commands.add("options(repos = c(CRAN = \"https://cloud.r-project.org\"))");
    commands.add("setwd(\"" + ext.verifyDirFormat(prefDir.getAbsolutePath()) + "\")");
    commands.add("if (!require(markdown)) {");
    commands.add("  install.packages(\"markdown\")");
    commands.add("}");
    commands.add("if (!require(MendelianRandomization)) {");
    commands.add("  install.packages(\"MendelianRandomization\")");
    commands.add("}");
    commands.add("library(MendelianRandomization)");
    commands.add("data <- read.table(\"" + MARKER_REGRESSION_PREFIX + pheno
                 + MARKER_REGRESSION_EXTEN + "\", sep=\"\\t\", header=T, stringsAsFactors=F)");
    commands.add("filtered <- data[data$BETA != \".\" & data$SE != \"0\",]");
    commands.add("bx <- filtered$META_BETA");
    commands.add("bxse <- filtered$META_SE");
    commands.add("by <- as.numeric(filtered$BETA)");
    commands.add("byse <- as.numeric(filtered$SE)");
    commands.add("input <- mr_input(bx=bx, bxse=bxse, by=by, byse=byse, snps=filtered$MARKERNAME, exposure=\"DepVar\", outcome=\"OutVar\")");
    commands.add("results <- mr_allmethods(input)");
    commands.add("print(results)");
    commands.add("write.table(results$Values, \"" + MENDELIAN_RANDOMIZATION_PREFIX + pheno
                 + MENDELIAN_RANDOMIZATION_EXTEN + "\", sep=\",\", row.names=F, quote=F)");
    Files.writeIterable(commands, ext.verifyDirFormat(prefDir.getAbsolutePath())
                                  + MENDELIAN_RANDOMIZATION_PREFIX + pheno + ".R");

    String setEnvCmd = Files.isWindows() ? "set" : "export";
    commands = new ArrayList<>();
    commands.add("cd " + prefDir.getAbsolutePath());
    if (rLibsDir != null) {
      commands.add(setEnvCmd + " R_LIBS=" + rLibsDir);
      commands.add(setEnvCmd + " R_LIBS_SITE=" + rLibsDir);
    } else {
      log.report("No R library directory specified, MendelianRandomization library will be installed if not present in default R library directory.");
    }
    commands.add("Rscript " + MENDELIAN_RANDOMIZATION_PREFIX + pheno + ".R");

    String runScript = ext.verifyDirFormat(prefDir.getAbsolutePath())
                       + MENDELIAN_RANDOMIZATION_PREFIX + pheno
                       + (Files.isWindows() ? ".bat" : ".sh");
    Files.writeIterable(commands, runScript);
    Files.chmod(runScript);
    return runScript;
  }

  private String getDirPath(Study study, MetaFile metaFile, Constraint constr) {
    return study.studyDir + metaFile.metaRoot + "/" + constr.analysisString + "/";
  }

  public static final String COMMAND_GENESCORE = "geneScore";

  public static void fromParameters(String filename, Logger log) {
    String dir = ext.verifyDirFormat(new File(ext.parseDirectoryOfFile(filename)).getAbsolutePath());
    List<String> params;
    params = Files.parseControlFile(filename, COMMAND_GENESCORE,
                                    new String[] {ARG_WORKING_DIR + dir,
                                                  ARG_INDEX_THRESH + DEFAULT_INDEX_THRESHOLD,
                                                  ARG_WINDOW_SIZE + DEFAULT_WINDOW_MIN_SIZE_PER_SIDE,
                                                  ARG_WINDOW_EXT + DEFAULT_WINDOW_EXTENSION_THRESHOLD},
                                    log);
    if (params == null) {
      setupCRF(filename, log);
      return;
    }
    if (params != null) {
      params.add("log=" + log.getFilename());
      main(ArrayUtils.toStringArray(params));
    }
  }

  private static void setupCRF(String crfFile, Logger log) {
    String dir = ext.parseDirectoryOfFile(crfFile);
    String[] files = (new File(dir)).list();
    HashSet<String> potentialRoots = new HashSet<>();
    HashSet<String> preprocess = new HashSet<>();
    HashSet<String> metas = new HashSet<>();
    HashSet<String> phenos = new HashSet<>();
    for (String f : files) {
      if (f.endsWith(".bed") || f.endsWith(".bim") || f.endsWith(".fam")) {
        potentialRoots.add(ext.rootOf(f));
      } else if (f.endsWith(".result") || f.endsWith(".results")) {
        preprocess.add(f);
      } else if (f.endsWith(".meta")) {
        metas.add(f);
      } else if (f.endsWith(".pheno")) {
        phenos.add(f);
      }
    }

    HashSet<String> validData = new HashSet<>();
    for (String s : potentialRoots) {
      if (PSF.Plink.allFilesExist(dir + s, true)) {
        validData.add(s);
      }
    }

    createDataFile(dir, validData, phenos);

    log.report("Processing " + preprocess.size() + " results files...");
    preprocessDataFiles(preprocess.toArray(new String[preprocess.size()]), POPULATION.ALL,
                        GenomeBuild.HG19, GenomeBuild.HG19, log);
    log.report("Done!");
  }

  private static void createDataFile(String baseDir, HashSet<String> plinkRoots,
                                     HashSet<String> phenoFiles) {
    String dir = baseDir + "plink/";
    new File(dir).mkdir();
    PrintWriter writer = Files.getAppropriateWriter(dir + "data.txt");
    for (String plinkRoot : plinkRoots) {
      writer.println(plinkRoot + "\t" + getFull(baseDir + PSF.Plink.getBED(plinkRoot)) + "\t"
                     + getFull(baseDir + PSF.Plink.getBIM(plinkRoot)) + "\t"
                     + getFull(baseDir + PSF.Plink.getFAM(plinkRoot)));
    }
    writer.flush();
    writer.close();

    for (String s : phenoFiles) {
      Files.copyFile(baseDir + s, dir + s);
    }
  }

  private static String getFull(String file) {
    return (new File(file)).getAbsolutePath();
  }

  public static void main(String[] args) {

    String workDir = null;
    String logFile = "GeneScorePipeline.log";

    float[] iT = new float[] {DEFAULT_INDEX_THRESHOLD};
    int[] mZ = new int[] {DEFAULT_WINDOW_MIN_SIZE_PER_SIDE};
    float[] wT = new float[] {DEFAULT_WINDOW_EXTENSION_THRESHOLD};
    double mT = DEFAULT_MIN_MISS_THRESH;
    int cmac = 5;

    String[] processList = null;
    boolean process = false;
    boolean runMetaHW = true;
    POPULATION pop = POPULATION.ALL;
    GenomeBuild build = GenomeBuild.HG19;
    GenomeBuild toBuild = GenomeBuild.HG19;
    String rLibsDir = null;
    boolean plotOddsRatio = false;
    boolean overwrite = false;

    String usage = "\n"
                   + "GeneScorePipeline is a convention-driven submodule.  It relies on a standard folder structure and file naming scheme:\n"
                   + "\tThe directory and file structure must conform to the following:\n"
                   + "\t\t>Root Directory ['workDir' argument]\n" + "\t\t\t>SNP Effect files:\n"
                   + "\t\t\t\t-Effect files must end with '.meta'.\n"
                   + "\t\t\t\t-Effect files may be hand-constructed, or may be generated with the 'preprocess' command from a .xln file\n"
                   + "\t\t\t\t-Effect files contain, at minimum, SNP, Freq, P-value, and Beta/Effect, and, if created with the preprocessor, will include any additional information present in the .xln file\n"
                   + "\t\t\t\t-HitWindows analysis will be run on SNP Effect files, with results being used in regression analysis; the\n"
                   + "\t\t\t\t\tadditional arguments to GeneScorePipeline affect only the HitWindows processing.\n"
                   + "\t\t\t>Covariate files:\n"
                   + "\t\t\t\t-Covariate files must end with '.covar'.\n"
                   + "\t\t\t\t-Covariate data will be added to ALL analyses.\n"
                   + "\t\t\t\t-Covariate files contain, at minimum, two ID columns (FID and IID), and any number of data columns.\n"
                   + "\t\t\t>Data Source Directory 1\n"
                   + "\t\t\t\t>data.txt file [defines location of data, which may be in an arbitrary location in the filesystem]\n"
                   + "\t\t\t\t\tExample:\n"
                   + "\t\t\t\t\tdataLabel1\tfullPathDataFile1\tFullPathMapFile1\tFullPathIdFile1\n"
                   + "\t\t\t\t\tdataLabel2\tfullPathDataFile2\tFullPathMapFile2\tFullPathIdFile2\n"
                   + "\t\t\t\t\tdataLabel3\tdir1\tdataFileExt1\tmapFileExt1\tidFile3\n"
                   + "\t\t\t\t\tdataLabel4\tdir2\tdataFileExt2\tmapFileExt2\tidFile4\n"
                   + "\t\t\t\t\t>For data files that contain map or ID info, 'FullPathMapFile' / 'FullPathIdFile' / 'mapFileExt1' / 'idFile' can be replaced with a period ('.').\n"
                   + "\t\t\t\t>Phenotype files\n"
                   + "\t\t\t\t\t-Phenotype files must end with '.pheno'.\n"
                   + "\t\t\t\t\t-All '.pheno' files will be included as a separate regression analysis\n"
                   + "\t\t\t\t\t-Phenotype files contain, at minimum, two ID columns (FID and IID), a dependent variable column, and any number of covariate columns.\n"
                   + "\t\t\t\t\t\t[Note: if data is in PLINK format and contains valid phenotype status information, an AFFECTED.PHENO file will be created]\n"
                   + "\t\t\t>Data Source Directory 2\n" + "\t\t\t\t>data.txt file\n"
                   + "\t\t\t\t>Pheno3.pheno file\n" + "\t\t\t>...\n" + "\n" + "\n"
                   + "gwas.GeneScorePipeline requires 1+ arguments\n"
                   + "   (1) Pre-process data files (i.e. process=path/to/file1.xln,path/to/file2.xln (not the default)) \n"
                   + "  OR\n"
                   + "   (1) Metastudy directory root, containing subdirectories for each study (i.e. workDir=C:/ (not the default))\n"
                   + "       OPTIONAL:\n"
                   + "   (2) p-value threshold (or comma-delimited list) for index SNPs (i.e. "
                   + ARG_INDEX_THRESH + DEFAULT_INDEX_THRESHOLD + " (default))\n"
                   + "   (3) minimum num bp per side of window (or comma delimited list) (i.e. "
                   + ARG_WINDOW_SIZE + DEFAULT_WINDOW_MIN_SIZE_PER_SIDE + " (default))\n"
                   + "   (4) p-value threshold to extend the window (or comma delimited list) (i.e. "
                   + ARG_WINDOW_EXT + DEFAULT_WINDOW_EXTENSION_THRESHOLD + " (default))\n"
                   + "   (5) minimum ratio of missing data for an individual's gene loading score to be included in the final analysis (i.e. "
                   + ARG_MISS_THRESH + DEFAULT_MIN_MISS_THRESH + " (default))\n"
                   + "   (6) Enable/Disable HitWindows on input '.meta' files (only applicable if #snps > 1000) (i.e. runMetaHW="
                   + runMetaHW + " (default))\n"

                   + "   (7) Directory of R libraries to generate Mendelian Randomization script from .meta SNPs (i.e. "
                   + ARG_R_LIBS_DIR + "/panfs/roc/msisoft/R/3.3.3/ (not activated by default))\n"

                   + "   (8) Flag to plot Odds Ratios isntead of Beta in Forest Plots (i.e. "
                   + FLAG_PLOT_ODDS_RATIO + " (not activated by default))\n"
                   + "   (9) Genome Build (i.e. build=" + build.name() + " (default))\n"
                   + "  (10) If pre-processing (with process= argument), request a LiftOver to a different genome build (liftBuild="
                   + toBuild.name() + " (default))\n" +

                   // " (8) Number of threads to use for computation (i.e. threads=" + threads + "
                   // (default))\n" +
                   "";

    boolean fail = false;
    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith(ARG_WORKING_DIR)) {
        workDir = arg.split("=")[1];
      } else if (arg.startsWith("process=")) {
        processList = arg.split("=")[1].split(",");
        process = true;
      } else if (arg.startsWith("pop=")) {
        pop = POPULATION.valueOf(arg.split("=")[1]);
      } else if (arg.startsWith("build=")) {
        build = GenomeBuild.valueOf(arg.split("=")[1].toUpperCase());
      } else if (arg.startsWith("liftBuild=")) {
        toBuild = GenomeBuild.valueOf(arg.split("=")[1].toUpperCase());
      } else if (arg.startsWith(ARG_INDEX_THRESH)) {
        String[] lst = arg.split("=")[1].split(",");
        int cntValid = 0;
        for (String poss : lst) {
          if (ext.isValidDouble(poss)) {
            cntValid++;
          }
        }
        iT = new float[cntValid];
        int ind = 0;
        for (String poss : lst) {
          if (ext.isValidDouble(poss)) {
            iT[ind] = Float.parseFloat(poss);
            ind++;
          }
        }
      } else if (arg.startsWith(ARG_WINDOW_SIZE)) {
        String[] lst = arg.split("=")[1].split(",");
        int cntValid = 0;
        for (String poss : lst) {
          if (ext.isValidDouble(poss)) {
            cntValid++;
          }
        }
        mZ = new int[cntValid];
        int ind = 0;
        for (String poss : lst) {
          if (ext.isValidInteger(poss)) {
            mZ[ind] = Integer.parseInt(poss);
            ind++;
          }
        }
      } else if (arg.startsWith(ARG_WINDOW_EXT)) {
        String[] lst = arg.split("=")[1].split(",");
        int cntValid = 0;
        for (String poss : lst) {
          if (ext.isValidDouble(poss)) {
            cntValid++;
          }
        }
        wT = new float[cntValid];
        int ind = 0;
        for (String poss : lst) {
          if (ext.isValidDouble(poss)) {
            wT[ind] = Float.parseFloat(poss);
            ind++;
          }
        }
      } else if (arg.startsWith(ARG_MISS_THRESH)) {
        mT = ext.parseDoubleArg(arg);
      } else if (arg.startsWith("runMetaHW=")) {
        runMetaHW = ext.parseBooleanArg(arg);
      } else if (arg.startsWith(ARG_R_LIBS_DIR)) {
        rLibsDir = ext.parseStringArg(arg);
      } else if (arg.equals(FLAG_PLOT_ODDS_RATIO)) {
        plotOddsRatio = true;
      } else if (arg.startsWith("log=")) {
        logFile = arg.split("=")[1];
      } else if (arg.equals("-overwrite")) {
        overwrite = true;
      } else if (arg.startsWith("cmac=")) {
        cmac = ext.parseIntArg(arg);
      } else {
        fail = true;
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (fail || args.length == 0) {
      System.err.println(usage);
      System.exit(1);
    }

    Logger log = new Logger(logFile, true);

    if (process) {
      preprocessDataFiles(processList, pop, build, toBuild, log);
      return;
    }

    File dir = new File(workDir);
    if (!dir.isDirectory()) {
      System.err.println("Error - argument 'workDir' must be a valid directory");
      System.exit(1);
    }
    GeneScorePipeline gsp = new GeneScorePipeline(workDir, iT, mZ, wT, mT, cmac, runMetaHW,
                                                  rLibsDir, plotOddsRatio, build, overwrite, log);
    gsp.runPipeline();
  }

}
