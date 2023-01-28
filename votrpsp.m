
clear all;
load votdata.txt;

d = votdata(:,1);
d2 = votdata(:,2);
d3 = votdata(:,3);
hinc = votdata(:,4);
exp = votdata(:,5);
free = ones(length(exp),1) - exp;

load beta_revised_het.txt;
b = beta_revised_het;
b0 = [b(1:6); b(28); b(33)];
load hess_revised_het.txt;
v = inv(hess_revised_het);
v0 = [v(:,1:6), v(:,28), v(:,33)];
v0 = [v0(1:6,:); v0(28,:); v0(33,:)];

cons = ones(length(d),1);

vot_50 = zeros(5000,1);
vot_het = zeros(5000,1);
vor_50 = zeros(5000,1);
vor_het = zeros(5000,1);

vote_50 = zeros(5000,1);
vote_het = zeros(5000,1);
vore_50 = zeros(5000,1);
vore_het = zeros(5000,1);

votf_50 = zeros(5000,1);
votf_het = zeros(5000,1);
vorf_50 = zeros(5000,1);
vorf_het = zeros(5000,1);

for i = 1:5000;
	bt = b0 + chol(v0)' * randn(length(b0),1);
	denom = bt(1) * cons + bt(2) * hinc;
	num1 = bt(3) * d + bt(4) * d2 + bt(5) * d3 + bt(7)*randn(length(d),1);
	num2 = bt(6) * cons + bt(6) * bt(8) * randn(length(d),1); 
	vot = num1*60 ./ denom; 
	vor = num2*60 ./ denom;
	
	vote = nonzeros(vot .* exp);
	votf = nonzeros(vot .* free); 

	vore = nonzeros(vor .* exp);
	vorf = nonzeros(vor .* free); 
	
	vot_50(i,1) = prctile(vot, 50);
	vot_het(i,1) = prctile(vot, 75) - prctile(vot, 25);
	vor_50(i,1) = prctile(vor, 50);
	vor_het(i,1) = prctile(vor,75) - prctile(vor, 25);

	vote_50(i,1) = prctile(vote, 50);
	vote_het(i,1) = prctile(vote, 75) - prctile(vote, 25);
	vore_50(i,1) = prctile(vore, 50);
	vore_het(i,1) = prctile(vore,75) - prctile(vore, 25);

	votf_50(i,1) = prctile(votf, 50);
	votf_het(i,1) = prctile(votf, 75) - prctile(votf, 25);
	vorf_50(i,1) = prctile(vorf, 50);
	vorf_het(i,1) = prctile(vorf,75) - prctile(vorf, 25);
end;


