import React, { useState } from 'react';

function App() {
    const [uniprotId, setUniprotId] = useState('');
    const [smiles, setSmiles] = useState('');
    const [loading, setLoading] = useState(false);
    const [result, setResult] = useState(null);
    const [error, setError] = useState(null);

    const handleSubmit = async (e) => {
        e.preventDefault();
        setLoading(true);
        setError(null);
        setResult(null);

        try {
            const response = await fetch('http://localhost:5000/api/dock', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify({ uniprot_id: uniprotId, smiles: smiles }),
            });

            if (!response.ok) {
                throw new Error(`Error: ${response.statusText}`);
            }

            const data = await response.json();
            setResult(data);
        } catch (err) {
            setError(err.message);
        } finally {
            setLoading(false);
        }
    };

    return (
        <div className="min-h-screen bg-gray-50 flex flex-col items-center py-12 px-4 sm:px-6 lg:px-8 font-sans">
            <div className="max-w-3xl w-full text-center mb-12">
                <h1 className="text-4xl font-extrabold text-gray-900 tracking-tight sm:text-5xl">
                    Enzyme-Substrate Simulator
                </h1>
                <p className="mt-4 text-lg text-gray-500">
                    Swiss-style precision docking on your local machine.
                </p>
            </div>

            <div className="max-w-md w-full bg-white p-8 border border-gray-200 shadow-sm rounded-lg">
                <form onSubmit={handleSubmit} className="space-y-6">
                    <div>
                        <label htmlFor="uniprot" className="block text-sm font-medium text-gray-700">
                            UniProt ID
                        </label>
                        <input
                            type="text"
                            id="uniprot"
                            className="mt-1 block w-full px-3 py-2 border border-gray-300 rounded-md shadow-sm focus:outline-none focus:ring-primary focus:border-primary sm:text-sm"
                            placeholder="e.g. P00533"
                            value={uniprotId}
                            onChange={(e) => setUniprotId(e.target.value)}
                            required
                        />
                    </div>

                    <div>
                        <label htmlFor="smiles" className="block text-sm font-medium text-gray-700">
                            SMILES String
                        </label>
                        <textarea
                            id="smiles"
                            rows={3}
                            className="mt-1 block w-full px-3 py-2 border border-gray-300 rounded-md shadow-sm focus:outline-none focus:ring-primary focus:border-primary sm:text-sm"
                            placeholder="e.g. CC(=O)Oc1ccccc1C(=O)O"
                            value={smiles}
                            onChange={(e) => setSmiles(e.target.value)}
                            required
                        />
                    </div>

                    <button
                        type="submit"
                        disabled={loading}
                        className={`w-full flex justify-center py-2 px-4 border border-transparent rounded-md shadow-sm text-sm font-medium text-white bg-primary hover:bg-green-800 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-primary ${loading ? 'opacity-75 cursor-not-allowed' : ''
                            }`}
                        style={{ backgroundColor: '#065F46' }} // Explicit inline for safety if tailwind fails
                    >
                        {loading ? 'Processing...' : 'Run Simulation'}
                    </button>
                </form>

                {error && (
                    <div className="mt-6 p-4 bg-red-50 border border-red-200 rounded-md">
                        <p className="text-sm text-red-600">{error}</p>
                    </div>
                )}
            </div>

            {result && (
                <div className="max-w-3xl w-full mt-12 grid grid-cols-1 gap-6 sm:grid-cols-2">
                    <div className="bg-white p-6 border border-gray-200 shadow-sm rounded-lg">
                        <h3 className="text-lg font-medium text-gray-900 border-b border-gray-100 pb-2 mb-4">Docking Affinity</h3>
                        <div className="text-5xl font-bold text-gray-900">
                            {result.affinity} <span className="text-xl font-normal text-gray-500">kcal/mol</span>
                        </div>
                    </div>

                    <div className="bg-white p-6 border border-gray-200 shadow-sm rounded-lg">
                        <h3 className="text-lg font-medium text-gray-900 border-b border-gray-100 pb-2 mb-4">Hydrogen Bonds</h3>
                        <div className="text-5xl font-bold text-primary" style={{ color: '#065F46' }}>
                            {result.h_bonds}
                        </div>
                        <p className="text-sm text-gray-500 mt-2">Interactions detected (&#60; 3.5&#197;)</p>
                    </div>

                    <div className="bg-white p-6 border border-gray-200 shadow-sm rounded-lg sm:col-span-2">
                        <h3 className="text-lg font-medium text-gray-900 border-b border-gray-100 pb-2 mb-4">Interaction Details</h3>
                        {result.h_bond_details.length > 0 ? (
                            <div className="overflow-x-auto">
                                <table className="min-w-full divide-y divide-gray-200 text-sm">
                                    <thead>
                                        <tr>
                                            <th className="px-4 py-2 text-left font-medium text-gray-500">Ligand Atom</th>
                                            <th className="px-4 py-2 text-left font-medium text-gray-500">Receptor Residue</th>
                                            <th className="px-4 py-2 text-left font-medium text-gray-500">Distance (&#197;)</th>
                                        </tr>
                                    </thead>
                                    <tbody className="divide-y divide-gray-200">
                                        {result.h_bond_details.map((bond, idx) => (
                                            <tr key={idx}>
                                                <td className="px-4 py-2 text-gray-900">{bond.ligand_atom}</td>
                                                <td className="px-4 py-2 text-gray-900">{bond.residue} ({bond.receptor_atom})</td>
                                                <td className="px-4 py-2 text-gray-900">{bond.distance.toFixed(2)}</td>
                                            </tr>
                                        ))}
                                    </tbody>
                                </table>
                            </div>
                        ) : (
                            <p className="text-gray-500 italic">No hydrogen bonds detected.</p>
                        )}
                    </div>

                    <div className="bg-gray-50 p-4 border border-gray-200 rounded text-xs font-mono text-gray-600 sm:col-span-2 overflow-auto max-h-64">
                        <pre>{result.stdout}</pre>
                    </div>
                </div>
            )}
        </div>
    );
}

export default App;
