/**
 * This worker contains all instructions to run the different steps that are required by the Peptonizer. All of these
 * functions are implemented in the same Worker instance and are loaded at the same time into memory instead of into
 * smaller separate workers. Since all of these workers require a running Pyodide instance, it's more efficient
 * (memory-wise) to bundle all Python functionality into one bigger worker.
 *
 * @author Pieter Verschaffelt
 */

import {loadPyodide, PyodideInterface} from 'pyodide';
import peptonizerWhlBase64 from "./lib/peptonizer-0.1-py3-none-any.base64.whl?raw";
import {
    ClusterTaxaTaskData,
    ClusterTaxaTaskDataResult, ComputeGoodnessDataResult, ComputeGoodnessTaskData,
    ExecutePepgmTaskData,
    ExecutePepgmTaskDataResult,
    FetchUnipeptTaxonTaskData,
    FetchUnipeptTaxonTaskResult,
    GenerateGraphTaskData,
    GenerateGraphTaskDataResult,
    InputEventData,
    OutputEventData,
    PerformTaxaWeighingTaskData,
    PerformTaxaWeighingTaskResult,
    ResultType,
    WorkerTask
} from "./PeptonizerWorkerTypes.ts";
import init, { 
    perform_taxa_weighing_wasm, 
    execute_pepgm_wasm, 
    fetch_unipept_taxa_wasm, 
    generate_pepgm_graph_wasm, 
    cluster_taxa_wasm,
    compute_goodness_wasm
} from "../../pkg/peptonizer_rust.js";

import fetchUnipeptTaxonPythonCode from "./lib/fetch_unipept_taxon_info.py?raw";
import performTaxaWeighingPythonCode from "./lib/perform_taxa_weighing.py?raw";
import generateGraphPythonCode from "./lib/generate_pepgm_graph.py?raw";
import executePepgmPythonCode from "./lib/execute_pepgm.py?raw";
import clusterTaxaPythonCode from "./lib/cluster_taxa.py?raw";
import computeGoodnessPythonCode from "./lib/compute_goodness.py?raw";
import { log, timeout } from 'async';

interface DedicatedWorkerGlobalScope {
    pyodide: PyodideInterface;
    postMessage: (message: OutputEventData) => void;
    submitPepgmProgress: (progressType: "graph" | "max_residual" | "iteration", currentValue: number, maxValue: number, workerId: number) => void;
}

declare const self: DedicatedWorkerGlobalScope & typeof globalThis;

async function loadPyodideAndPackages(): Promise<void> {
    self.pyodide = await loadPyodide({
        indexURL: 'https://cdn.jsdelivr.net/pyodide/v0.26.3/full/'
    });

    // Load all packages into the Pyodide runtime environment that are required by the Peptonizer
    await self.pyodide.loadPackage([
        'numpy',
        'scipy',
        'networkx',
        'pandas',
        'micropip',
        'requests',
        'openssl'
    ]);
    // Use the imported .whl file URL directly with micropip
    await self.pyodide.runPythonAsync(`
        import base64
        from pathlib import Path
        
        import micropip

        await micropip.install('rbo')

        # Decode base64 string to binary and write to a temporary file
        wheel_data = "${peptonizerWhlBase64.replace(/\s+/g, '')}"
        wheel_binary = base64.b64decode(wheel_data)

        # Define a temporary path for the .whl file
        wheel_path = Path("/tmp/peptonizer-0.1-py3-none-any.whl")
        wheel_path.write_bytes(wheel_binary)

        # Install the wheel package
        await micropip.install("emfs:///tmp/peptonizer-0.1-py3-none-any.whl")

        # Clean up by deleting the temporary file
        wheel_path.unlink()
    `);
}

async function fetchUnipeptTaxonInformation(data: FetchUnipeptTaxonTaskData): Promise<FetchUnipeptTaxonTaskResult> {
    console.time("Execution time fetching Unipept information");
    // Set inputs for the Python code
    /*self.pyodide.globals.set('peptides_scores', data.peptidesScores);
    self.pyodide.globals.set('rank', data.rank);
    self.pyodide.globals.set('taxon_query', data.taxonQuery);

    const unipeptJson = await self.pyodide.runPythonAsync(fetchUnipeptTaxonPythonCode);*/
    
    let score_keys = [...data.peptidesScores.keys()];
    let peptidesScores = JSON.stringify(score_keys);
    let taxonQuery = JSON.stringify(data.taxonQuery);

    const unipeptJson = fetch_unipept_taxa_wasm(peptidesScores, data.rank, taxonQuery);
    
    console.timeEnd("Execution time fetching Unipept information");

    return { unipeptJson };
}

async function performTaxaWeighing(data: PerformTaxaWeighingTaskData): Promise<PerformTaxaWeighingTaskResult> {
    console.time("Execution time taxa weiging");
    
    let peptidesTaxa = JSON.stringify(Object.fromEntries(data.peptidesTaxa));
    let peptidesScores = JSON.stringify(Object.fromEntries(data.peptidesScores));
    let peptidesCounts = JSON.stringify(Object.fromEntries(data.peptidesCounts));

    const [sequenceScoresCsv, taxaWeightsCsv] = perform_taxa_weighing_wasm(peptidesTaxa, peptidesScores, peptidesCounts, data.taxaInGraph, "species");
    
    // Set inputs for the Python code
    /*self.pyodide.globals.set('peptides_taxa', data.peptidesTaxa);
    self.pyodide.globals.set('peptides_scores', data.peptidesScores);
    self.pyodide.globals.set('peptides_counts', data.peptidesCounts);
    self.pyodide.globals.set('rank', data.rank);
    self.pyodide.globals.set('taxa_in_graph', data.taxaInGraph);

    // Fetch the Python code and execute it with Pyodide
    const [sequenceScoresCsv_py, taxaWeightsCsv_py] = await self.pyodide.runPythonAsync(performTaxaWeighingPythonCode);
    */

    console.timeEnd("Execution time taxa weiging");
    return {
        sequenceScoresCsv,
        taxaWeightsCsv
    };

    
}

async function generateGraph(data: GenerateGraphTaskData): Promise<GenerateGraphTaskDataResult> {
    console.time("Execution time generating graph");

    const graphXml = generate_pepgm_graph_wasm(data.taxaWeightsCsv);

    /*self.pyodide.globals.set('taxa_weights_csv', data.taxaWeightsCsv);
    const graphXml = await self.pyodide.runPythonAsync(generateGraphPythonCode);*/
    
    console.timeEnd("Execution time generating graph");

    return {
        graphXml
    };
}


async function executePepgm(data: ExecutePepgmTaskData, workerId: number): Promise<ExecutePepgmTaskDataResult> {
    console.time("Execution time PepGM");

    const taxonScoresJson = execute_pepgm_wasm(data.graphXml, data.alpha, data.beta, true, data.prior);

    /*self.pyodide.globals.set('graph', data.graphXml);
    self.pyodide.globals.set('alpha', data.alpha);
    self.pyodide.globals.set('beta', data.beta);
    self.pyodide.globals.set('prior', data.prior);
    self.pyodide.globals.set('worker_id', workerId);

    const taxonScoresJson = await self.pyodide.runPythonAsync(executePepgmPythonCode);
    */

    console.timeEnd("Execution time PepGM");
    return {
        taxonScoresJson
    };
}

async function clusterTaxa(data: ClusterTaxaTaskData): Promise<ClusterTaxaTaskDataResult> {
    console.time("Execution time clustering taxa");

    const clusteredTaxaWeightsCsv = cluster_taxa_wasm(data.graphXml, data.taxaWeightsCsv, data.similarityThreshold)

    /*self.pyodide.globals.set('graph', data.graphXml);
    self.pyodide.globals.set('taxa_weights_csv', data.taxaWeightsCsv);
    self.pyodide.globals.set('similarity_threshold', data.similarityThreshold);

    const clusteredTaxaWeightsCsv = await self.pyodide.runPythonAsync(clusterTaxaPythonCode);*/

    console.timeEnd("Execution time clustering taxa");
    return {
        clusteredTaxaWeightsCsv
    };
}

async function computeGoodness(data: ComputeGoodnessTaskData): Promise<ComputeGoodnessDataResult> {
    console.time("Execution time computing goodness");
    
    let peptonizerResults = JSON.stringify(Object.fromEntries(data.peptonizerResults));
    const goodness = compute_goodness_wasm(data.clusteredTaxaWeightsCsv, peptonizerResults);
    
    /*self.pyodide.globals.set('clustered_taxa_weights_csv', data.clusteredTaxaWeightsCsv);
    self.pyodide.globals.set('peptonizer_results', data.peptonizerResults);

    const goodness = await self.pyodide.runPythonAsync(computeGoodnessPythonCode);
    */
    console.log(goodness);

    console.timeEnd("Execution time computing goodness");
    return {
        goodness
    }
}

self.submitPepgmProgress = function(
    progressType: "graph" | "max_residual" | "iteration",
    currentValue: number,
    maxValue: number,
    workerId: number
) {
    const resultMessage: OutputEventData = {
        resultType: ResultType.PROGRESS,
        task: WorkerTask.EXECUTE_PEPGM,
        workerId: workerId,
        progressUpdate: {
            progressType,
            currentValue,
            maxValue
        }
    }

    self.postMessage(resultMessage);
}

// let pyodideReadyPromise: Promise<void> = loadPyodideAndPackages();

self.onmessage = async (event: MessageEvent<InputEventData>): Promise<void> => {
    try {
        // Make sure loading is done
        // await pyodideReadyPromise;
        await init();

        // Destructure the data from the event
        const eventData = event.data;

        let output: FetchUnipeptTaxonTaskResult | PerformTaxaWeighingTaskResult | GenerateGraphTaskDataResult | ExecutePepgmTaskDataResult | ClusterTaxaTaskDataResult | ComputeGoodnessDataResult | undefined;

        if (eventData.task === WorkerTask.FETCH_UNIPEPT_TAXON) {
            output = await fetchUnipeptTaxonInformation(eventData.input);
        } else if (eventData.task === WorkerTask.PERFORM_TAXA_WEIGHING) {
            output = await performTaxaWeighing(eventData.input);
        } else if (eventData.task === WorkerTask.GENERATE_GRAPH) {
            output = await generateGraph(eventData.input);
        } else if (eventData.task === WorkerTask.EXECUTE_PEPGM) {
            output = await executePepgm(eventData.input, eventData.workerId);
        } else if (eventData.task === WorkerTask.CLUSTER_TAXA) {
            output = await clusterTaxa(eventData.input);
        } else if (eventData.task === WorkerTask.COMPUTE_GOODNESS) {
            output = await computeGoodness(eventData.input);
        } else {
            throw new Error("Unknown task type passed to worker!");
        }

        if (!output) {
            throw new Error("No valid output was generated by worker!");
        }

        self.postMessage({
            resultType: ResultType.SUCCESSFUL,
            workerId: eventData.workerId,
            task: eventData.task,
            output
        });
    } catch (error: any) {
        self.postMessage({
            resultType: ResultType.FAILED,
            workerId: event.data.workerId,
            error: error.toString()
        });
    }
};
